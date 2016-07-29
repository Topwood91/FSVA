#include	"radix_sort.cpp"
#include	"ssw_cpp.h"
#include	<algorithm>
#include	<iostream>
#include	<stdint.h>
#include	<cstdlib>
#include	<cstring>
#include	<fstream>
#include	<sstream>
#include	<vector>
#include	<zlib.h>
#include	<cmath>
#include	<omp.h>
#include    <string>
using	namespace	std;
const	unsigned	mod=(1UL<<29)-1;
unsigned	kmer=31;
uint64_t	mask;
unsigned	max_occ;
const	unsigned	max_indel=32;
const	unsigned	buffer=1000000;
const   unsigned    qkmp=6;
const   unsigned    pairdis=1000;
unsigned    max_votes;
unsigned    readLen;
unsigned	r1er2;
uint8_t	code[256];
char	symbol[6]="ACGTN",	rcsymbol[6]="TGCAN";
void	make_code(void){
	for(size_t	i=0;	i<256;	i++)	code[i]=4;
	code['A']=code['a']=0;	code['C']=code['c']=1;	code['G']=code['g']=2;	code['T']=code['t']=3;
}

struct	Read{
	string	name,	seq,	qua,	cigar;
	uint32_t	tid,	pos,	as;
	char	strand,	mapq;
};

struct	Region{
	uint32_t	beg;
	uint32_t	end;
	uint32_t	pos;
	uint16_t	cov;
	char	str;
	bool	operator()(Region	X,	Region	Y){	return	X.cov<Y.cov;	}
};

class	fuwa{
private:
	StripedSmithWaterman::Aligner aligner;
	StripedSmithWaterman::Filter filter;
	string	ref,rgtag;
	vector<string>	name;
	vector<uint32_t>	offset,	keyv,	posv;
	vector<Read>	read;
	vector<Read>    read2;
	void	map_block(size_t	N, char mode);
	void	map_read(Read	&R);
	void    map_read2(Read	&R1, Read  &R2);
	void	query(string	&Q,	vector<Region>	&R,vector<Region> &R2,	char	S,	size_t	&M,	size_t	&N,bool ispe);
public:
	bool	sam_header(void);
	bool	load(const	char	*F);
	bool	fastq(const	char	*F, const char *F2);
};

bool	fuwa::load(const	char	*F ){
	aligner.Clear();	aligner.ReBuild(1,4,6,1);
	FILE	*f=fopen(F,	"rb");
	if(f==NULL)	return	false;
	fseek(f,	0,	SEEK_END);
	ref.reserve(ftell(f)+1);
	fclose(f);

	ifstream	fi(F);
	if(!fi){	cerr<<"fail to open "<<F<<'\n';	return	false;	}
	cerr<<"loading references\n";
	string	buf,	temp;
	while(!fi.eof()){
		getline(fi,	buf);	if(!buf.size())	continue;
		if(buf[0]=='>'){
			istringstream	si(buf.c_str()+1);
			string	t;	si>>t;	name.push_back(t);	offset.push_back(ref.size());
		}
		else{
			for(size_t	i=0;	i<buf.size();	i++)	if(buf[i]>33)	ref.push_back(*(code+buf[i]));
		}
	}
	fi.close();
	offset.push_back(ref.size());
	cerr<<ref.size()<<endl;
	string	fn=F;	fn+=".hash";
	fi.open(fn.c_str(),	ios::binary);

	uint32_t	n;
	fi.read((char*)&n,	4);
	posv.resize(n);
	fi.read((char*)&posv[0],	posv.size()*4);
	keyv.resize(mod+1);
	fi.read((char*)&keyv[0],	keyv.size()*4);
	fi.close();
	return	true;
}

bool	fuwa::fastq(const	char	*F1, const char *F2){
	gzFile	in1=gzopen(F1,	"rt");
	if(in1==Z_NULL)	return	false;
	gzFile in2=gzopen(F2,"rt");
	if(strlen(F2)>0 && in2==Z_NULL)return false;
	char	name[65536],	seq[65536],	qua[65536],	temp[65536];
	read.resize(buffer);
	read2.resize(buffer);
	size_t	reads1=0;
	size_t  reads2=0;
	while(!gzeof(in1)){
		if(gzgets(in1,	name,	65536)==NULL)	break;
		if(gzgets(in1,	seq,	65536)==NULL)	break;
		if(gzgets(in1,	temp,	65536)==NULL)	break;
		if(gzgets(in1,	qua,	65536)==NULL)	break;
		Read	&r=read[reads1%buffer];
		size_t	len=strlen(seq)-1;
		istringstream	si(name+1);
		si>>r.name;
		r.seq.assign(seq,	len);
		r.qua.assign(qua,	len);
		if(strlen(F2)>0)
        {
            if(gzgets(in2,	name,	65536)==NULL)	break;
            if(gzgets(in2,	seq,	65536)==NULL)	break;
            if(gzgets(in2,	temp,	65536)==NULL)	break;
            if(gzgets(in2,	qua,	65536)==NULL)	break;
            Read	&r=read2[reads2%buffer];
            size_t	len=strlen(seq)-1;
            istringstream	si(name+1);
            si>>r.name;
            r.seq.assign(seq,	len);
            r.qua.assign(qua,	len);
            reads2++;
        }
		if(++reads1%buffer==0){	map_block(buffer,strlen(F2)>0?'p':'s');	if(reads1==buffer)	cerr<<F1<<'\t'<<F2<<'\t';	cerr<<'=';	}
	}
	if(reads1%buffer){	map_block(reads1%buffer,strlen(F2)>0?'p':'s');	cerr<<'=';	}
	cerr<<'\n';
	gzclose(in1);
	if(strlen(F2)>0)gzclose(in2);
	return	true;
}

void	fuwa::map_block(size_t	N, char mode){
	if(mode=='s')
    {
        #pragma omp parallel for
        for(size_t	i=0;	i<N;	i++)	map_read(read[i]);
        for(size_t	i=0;	i<N;	i++){
            ostringstream	so;
            printf("%s\t",	read[i].name.c_str());
            if(read[i].strand=='*')	printf("4\t*\t0\t0\t*\t*\t0\t0\t");
            else	printf("%d\t%s\t%u\t%d\t%s\t*\t0\t0\t",	read[i].strand=='+'?0:16,	name[read[i].tid].c_str(),	read[i].pos,	(int)read[i].mapq,	read[i].cigar.c_str());
            if(read[i].strand=='-'){
                size_t	len=read[i].seq.size();
                for(size_t	j=0;	j<len;	j++)	putchar(*(rcsymbol+*(code+read[i].seq[len-1-j])));
                putchar('\t');
                for(size_t	j=0;	j<len;	j++)	putchar(read[i].qua[len-1-j]);
            }
            else	printf("%s\t%s",	read[i].seq.c_str(),	read[i].qua.c_str());
            if(read[i].strand!='*')	printf("\tAS:i:%u",	read[i].as);
            putchar('\n');
        }
    }
    else {
        #pragma omp parallel for
        for(size_t i=0;i<N;i++) map_read2(read[i],read2[i]);
        for(size_t	i=0;	i<N;	i++){
            ostringstream	so;

            //read1
            so<<read[i].name<<'\t';
            uint16_t	flag=0x1;
            if(read[i].mapq==0)	flag|=0x4;
            if(read2[i].mapq==0) flag|=0x8;
            if(!(flag&0x4)&&!(flag&0x8))	flag|=0x2;
            if(read[i].mapq&&read[i].strand=='-')	flag|=0x10;
            if(read2[i].mapq&&read2[i].strand=='-')	flag|=0x20;
            flag|=0x40;
            so<<flag;
            so<<'\t'<<(read[i].mapq==0?"*":name[read[i].tid]);
            so<<'\t'<<(read[i].mapq==0?0:read[i].pos);
            so<<'\t'<<(read[i].mapq==0?0:(int)read[i].mapq)<<'\t';

            if(read[i].mapq==0)	so<<"*";
            else so<<read[i].cigar;
            if(read2[i].mapq)   so<<'\t'<<name[read2[i].tid]<<'\t'<<read2[i].pos<<"\t0\t";
            else	so<<"\t*\t0\t0\t";
            if(read[i].strand=='-'){
                size_t	len=read[i].seq.size();
                for(size_t	j=0;	j<len;	j++)	so<<*(rcsymbol+*(code+read[i].seq[len-1-j]));
                so<<'\t';
                for(size_t	j=0;	j<len;	j++)	so<<read[i].qua[len-1-j];
            }
            else	so<<read[i].seq<<'\t'<<read[i].qua;
           // so<<'\t'<<rgtag<<'\n';
            cout<<so.str()<<endl;

            //read2
            so.str("");
            so<<read2[i].name<<'\t';
            flag=0x1;
            if(read2[i].mapq==0)	flag|=0x4;
            if(read[i].mapq==0) flag|=0x8;
            if(!(flag&0x4)&&!(flag&0x8))	flag|=0x2;
            if(read2[i].mapq&&read2[i].strand=='-')	flag|=0x10;
            if(read[i].mapq&&read[i].strand=='-')	flag|=0x20;
            flag|=0x40;
            so<<flag;
            so<<'\t'<<(read2[i].mapq==0?"*":name[read2[i].tid]);
            so<<'\t'<<(read2[i].mapq==0?0:read2[i].pos);
            so<<'\t'<<(read2[i].mapq==0?0:(int)read2[i].mapq)<<'\t';
            if(read2[i].mapq==0)	so<<"*";
            else so<<read2[i].cigar;
            if(read[i].mapq)   so<<'\t'<<name[read[i].tid]<<'\t'<<read[i].pos<<"\t0\t";
            else	so<<"\t*\t0\t0\t";
            if(read2[i].strand=='-'){
                size_t	len=read2[i].seq.size();
                for(size_t	j=0;	j<len;	j++)	so<<*(rcsymbol+*(code+read2[i].seq[len-1-j]));
                so<<'\t';
                for(size_t	j=0;	j<len;	j++)	so<<read2[i].qua[len-1-j];
            }
            else	so<<read2[i].seq<<'\t'<<read2[i].qua;
          //  so<<'\t'<<rgtag<<'\n';
            cout<<so.str()<<endl;
        }

    }

}

void	fuwa::map_read(Read	&R){
	size_t	len=R.seq.size(),	best=0,	next=0;
	string	fwd,	rev;	fwd.resize(len);	rev.resize(len);
	vector<Region>	vr,vr2;	vr.reserve(1024);
	for(size_t	i=0;	i<len;	i++){
		uint8_t	c=*(code+R.seq[i]);
		fwd[i]=c;	rev[len-1-i]=c==4?4:3-c;
	}
	query(fwd,	vr,vr2,	'+',	best,	next,false);	query(rev,	vr,vr2,	'-',	best,	next,false);
	if(!vr.size()){	R.strand='*';	return;	}
	if(vr.size()>1)	next=best;	R.mapq=fminf((best-next)*6,	60);
    size_t	i=rand()%vr.size();

	string	&s=vr[i].str=='+'?fwd:rev;
	unsigned	mismatch=0;
	if(vr[i].pos+len<=ref.size()){
		const	char	*p=ref.c_str()+vr[i].pos;
		for(size_t	j=0;	j<len;	j++)	if(s[j]!=p[j])	mismatch++;
	}

	else	mismatch = 2;

	if(mismatch<2){
		R.pos=vr[i].pos;
		R.tid=0;	for(size_t	j=0;	j<name.size();	j++)	if(offset[j+1]>R.pos){	R.tid=j;	break;	}
		R.pos=R.pos-offset[R.tid]+1;	R.strand=vr[i].str;
		char	temp[1024];	sprintf(temp,	"%uM",	(unsigned)len);	R.cigar=temp;
		R.as=len-5*mismatch;
	}
	else{
		StripedSmithWaterman::Alignment a;
		aligner.Align(s.c_str(),	len,	ref.c_str()+vr[i].beg,	vr[i].end-vr[i].beg,	filter, &a);
		R.pos=vr[i].beg+a.ref_begin;
		R.tid=0;	for(size_t	j=0;	j<name.size();	j++)	if(offset[j+1]>R.pos){	R.tid=j;	break;	}
		R.pos=R.pos-offset[R.tid]+1;	R.strand=vr[i].str;	R.cigar=a.cigar_string;
		R.as=a.sw_score;
	}
}

void fuwa::map_read2(Read	&R1, Read  &R2) {
	size_t	len = R1.seq.size(), best = 0, next = 0;
	size_t  len2=R2.seq.size(),best2=0,next2=0;
	string	fwd, rev;	fwd.resize(len);	rev.resize(len);
	string  fwd2,rev2;  fwd2.resize(len2);  rev2.resize(len2);
	vector<Region>	vr, vr2;	vr.reserve(1024);   vr2.reserve(1024);
	vector<Region>  vrt,vrt2;   vrt.reserve(1024);  vrt2.reserve(1024);
	Region choosen1,choosen2;
	for (size_t i = 0; i < len; i++) {
		uint8_t	c = *(code + R1.seq[i]);
		fwd[i] = c;	rev[len - 1 - i] = c == 4 ? 4 : 3 - c;
	}
	for(size_t i=0;i<len2;i++){
        uint8_t c=*(code+R2.seq[i]);
        fwd2[i]=c;rev2[len2-1-i]=c==4?4:3-c;
	}
	query(fwd, vr, vr2, '+', best, next, true);	query(rev, vr, vr2, '-', best, next, true);
	query(fwd2,vrt,vrt2,'+',best2,next2,true);query(rev2,vrt,vrt2,'-',best2,next2,true);
	if (!vr.size()) { R1.strand = '*';	return; }
	if(!vrt.size()){R2.strand='*';return;}

	choosen1 = vr[rand() % vr.size()];
	choosen2=vrt[rand()%vrt.size()];
	if (vr.size()>1 && vrt.size()>1) {
		next = best;
		next2=best2;
	}
	else if(vr.size()==1 && vrt.size()>1){
        next2=best2;
        bool flag=false;
        for(size_t i=0;i<vrt.size();i++){
            if((choosen1.str!=vrt[i].str) && ((choosen1.pos>vrt[i].pos&&choosen1.pos-vrt[i].pos<=pairdis) || (choosen1.pos<=vrt[i].pos && vrt[i].pos-choosen1.pos<=pairdis))){
                next*=2;best*=2;
                next2=best2-1;
                choosen2=vrt[i];
                flag=true;
                break;
            }
        }
        if(!flag){
            for(size_t i=0;i<vrt2.size();i++){
                if((choosen1.str!=vrt2[i].str) && ((choosen1.pos>vrt2[i].pos&&choosen1.pos-vrt2[i].pos<=pairdis) || (choosen1.pos<=vrt2[i].pos && vrt2[i].pos-choosen1.pos<=pairdis))){
                    choosen2=vrt2[i];
                    break;
                }
            }
        }
	}
	else if(vr.size()>1 && vrt.size()==1){
        next=best;
        bool flag=false;
        for(size_t i=0;i<vr.size();i++){
            if((choosen2.str!=vr[i].str) && ((choosen2.pos>vr[i].pos&&choosen2.pos-vr[i].pos<=pairdis) || (choosen2.pos<=vr[i].pos && vr[i].pos-choosen2.pos<=pairdis))){
                next2*=2;best2*=2;
                next=best-1;
                choosen1=vr[i];
                flag=true;
                break;
            }
        }
        if(!flag){
            for(size_t i=0;i<vr2.size();i++){
                if((choosen2.str!=vr2[i].str) && ((choosen2.pos>vr2[i].pos&&choosen2.pos-vr2[i].pos<=pairdis) || (choosen2.pos<=vr2[i].pos && vr2[i].pos-choosen2.pos<=pairdis))){
                    choosen1=vr2[i];
                    flag=true;
                    break;
                }
            }
        }
	}
	else {
        if((choosen1.str!=choosen2.str) && ((choosen1.pos>choosen2.pos && choosen1.pos-choosen2.pos<=pairdis) || (choosen1.pos<=choosen2.pos && choosen2.pos-choosen1.pos<=pairdis))){
            best*=3;    best2*=3;
            next*=3;    next2*=3;
        }
	}
	R1.mapq = fminf((best - next) * 6, 60);
	R2.mapq = fminf((best2-next2)*6,60);
	string	&s = choosen1.str == '+' ? fwd : rev;
	string  &s2=choosen2.str=='+'?fwd2:rev2;
	unsigned	mismatch = 0;
	if (choosen1.pos + len <= ref.size()) {
		const	char	*p = ref.c_str() + choosen1.pos;
		for (size_t j = 0; j < len; j++)	if (s[j] != p[j])	mismatch++;
	}

	else	mismatch = 2;
	if (mismatch < 2) {
		R1.pos = choosen1.pos;
		R1.tid = 0;	for (size_t j = 0; j<name.size(); j++)	if (offset[j + 1]>R1.pos) { R1.tid = j;	break; }
		R1.pos = R1.pos - offset[R1.tid] + 1;	R1.strand = choosen1.str;
		char	temp[1024];	sprintf(temp, "%uM", (unsigned)len);	R1.cigar = temp;
		R1.as = len - 5 * mismatch;
	}

	else {
		StripedSmithWaterman::Alignment a;
		aligner.Align(s.c_str(), len, ref.c_str() + choosen1.beg, choosen1.end - choosen1.beg, filter, &a);
		R1.pos = choosen1.beg + a.ref_begin;
		R1.tid = 0;	for (size_t j = 0; j<name.size(); j++)	if (offset[j + 1]>R1.pos) { R1.tid = j;	break; }
		R1.pos = R1.pos - offset[R1.tid] + 1;	R1.strand = choosen1.str;	R1.cigar = a.cigar_string;
		R1.as = a.sw_score;
	}

	mismatch=0;
	if (choosen2.pos + len2 <= ref.size()) {
		const	char	*p = ref.c_str() + choosen2.pos;
		for (size_t j = 0; j < len2; j++)	if (s2[j] != p[j])	mismatch++;
	}

	else	mismatch = 2;
	if (mismatch < 2) {
		R2.pos = choosen2.pos;
		R2.tid = 0;	for (size_t j = 0; j<name.size(); j++)	if (offset[j + 1]>R2.pos) { R2.tid = j;	break; }
		R2.pos = R2.pos - offset[R2.tid] + 1;	R2.strand = choosen2.str;
		char	temp[1024];	sprintf(temp, "%uM", (unsigned)len2);	R2.cigar = temp;
		R2.as = len2 - 5 * mismatch;
	}

	else {
		StripedSmithWaterman::Alignment a;
		aligner.Align(s2.c_str(), len2, ref.c_str() + choosen2.beg, choosen2.end - choosen2.beg, filter, &a);
		R2.pos = choosen2.beg + a.ref_begin;
		R2.tid = 0;	for (size_t j = 0; j<name.size(); j++)	if (offset[j + 1]>R2.pos) { R2.tid = j;	break; }
		R2.pos = R2.pos - offset[R2.tid] + 1;	R2.strand = choosen2.str;	R2.cigar = a.cigar_string;
		R2.as = a.sw_score;
	}
}

void	fuwa::query(string	&Q,	vector<Region>	&R, vector<Region>	&R2, char	S,	size_t	&M,	size_t	&N, bool ispe){
	vector<uint32_t>	hit;	hit.reserve(1024);
	uint64_t	k=0;
	for(size_t	i=0;	i<kmer-1;	i++)	k=(k<<2)+(Q[i]==4?rand()&3:Q[i]);
	for(size_t	i=kmer-1;	i<Q.size();	i++){
		k=(k<<2)+(Q[i]==4?rand()&3:Q[i]);
		size_t	h=(k&mask)%mod,	b=keyv[h],	e=keyv[h+1];
		if(e-b>max_occ)	continue;
		for(size_t	j=b;	j<e;	j++)	hit.push_back(posv[j]+kmer>=1+i?posv[j]+kmer-1-i:0);
	}
	radix_sort(&hit[0],	&hit[0]+hit.size());
	for(size_t	i=0;	i<hit.size();){
		size_t	j;
		for(j=i+1;	j<hit.size()&&hit[j]==hit[i];	j++);
		Region	r;	r.cov=j-i;
		if(r.cov>=M&&r.cov>=2){
			if (r.cov > M) {N = M;	M = r.cov; R2=R;	R.clear(); }
			r.pos=hit[i];
			r.beg=r.pos>max_indel?r.pos-max_indel:0;
			r.end=r.pos+Q.size()+max_indel<ref.size()?r.pos+Q.size()+max_indel:ref.size();
			r.str=S;
			R.push_back(r);
		}
		else if(r.cov>=N && r.cov>=2){
            if(ispe){
                if(r.cov>N){N = r.cov; R2.clear();}
                r.pos=hit[i];
                r.beg=r.pos>max_indel?r.pos-max_indel:0;
                r.end=r.pos+Q.size()+max_indel<ref.size()?r.pos+Q.size()+max_indel:ref.size();
                r.str=S;
                R2.push_back(r);
            }
            else{
                if(r.cov>N)N=r.cov;
            }
		}
		i=j;
	}
}

bool	fuwa::sam_header(void){
	cout<<"@HD	VN:1.3	SO:coordinate\n";
	for(size_t	i=0;	i<name.size();	i++)	cout<<"@SQ\tSN:"<<name[i]<<'\t'<<"LN:"<<offset[i+1]-offset[i]<<'\n';
	cout<<"@PG	ID:fuwa	PN:fuwa	VN:0.0\n";
	return	true;
}

int	main(int	ac,	char	**av){
    max_occ=450;
    int threads=1;
    kmer=31;
	if(ac<3){	cerr<<"fsva [options] <ref.fa> <read1.fastq> [read2.fastq]\n";
        cerr<<"options:\n";
        cerr<<"\t-t INT number of threads. [1]\n";
        cerr<<"\t-u INT threshold of unrepresentative seed. [450]\n";
        cerr<<"\t-l INT length of seed. If you use this argument when making\n\t       index, you should use this argument here, and their\n\t       value should be equal. [31]\n";
        cerr<<"\t-r INT length of read. If you use this argument when making\n\t       index, you should use this argument here, and their\n\t       value should be equal. [150]\n\n";
	return	0;	}
	int opn=1;
	int kmer_temp=0;
	while(opn<ac){
        bool flag=false;
        if(av[opn][0]=='-'){
            if(av[opn][1]=='t'){threads=atoi(av[opn+1]);opn+=2;flag=true;}
            else if(av[opn][1]=='u'){max_occ=atoi(av[opn+1]);opn+=2;flag=true;}
            else if(av[opn][1]=='l'){kmer_temp=atoi(av[opn+1]);opn+=2;flag=true;}
            else if(av[opn][1]=='r'){readLen=atoi(av[opn+1]);opn+=2;flag=true;}
            else {
                cerr<<"fsva [options] <ref.fa> <read1.fastq> [read2.fastq]\n";
                cerr<<"options:\n";
                cerr<<"\t-t INT number of threads. [1]\n";
                cerr<<"\t-u INT threshold of unrepresentative seed. [450]\n";
                cerr<<"\t-l INT length of seed. If you use this argument when making\n\t       index, you should use this argument here, and their\n\t       value should be equal. [31]\n";
                cerr<<"\t-r INT length of read. If you use this argument when making\n\t       index, you should use this argument here, and their\n\t       value should be equal. [150]\n\n";
                return	0;
            }
        }
        if(!flag)break;
	}
	if(readLen>=175)kmer=32;
    else if(readLen>=150)kmer=31;
    else if(readLen>=125)kmer=30;
    else if(readLen>=75)kmer=21;
    else if(readLen>0)kmer=16;
    if(kmer_temp!=0)kmer=kmer_temp;
	mask=kmer==32?~0:(1ULL<<(kmer*2))-1;
	max_votes=(readLen-kmer+1)/8;
	omp_set_num_threads(threads);
	r1er2 = 0;
	make_code();
	fuwa	f;
	setvbuf(stdout,	NULL,	_IOFBF,	16*1024*1024);
	if(!f.load(av[opn]))	return	0;
	opn++;
	if(!f.sam_header())	return	0;
	size_t	t0=time(NULL);
	if(opn==ac-1) f.fastq(av[opn],"\0");
	else if(opn==ac-2)f.fastq(av[opn],av[opn+1]);
	else{
        cerr<<"fsva [options] <ref.fa> <read1.fastq> [read2.fastq]\n";
        cerr<<"options:\n";
        cerr<<"\t-t INT number of threads. [1]\n";
        cerr<<"\t-u INT threshold of unrepresentative seed. [450]\n";
        cerr<<"\t-l INT length of seed. If you use this argument when making\n\t       index, you should use this argument here, and their\n\t       value should be equal. [31]\n";
        cerr<<"\t-r INT length of read. If you use this argument when making\n\t       index, you should use this argument here, and their\n\t       value should be equal. [150]\n\n";
        return	0;
	}
	cerr << r1er2 << endl;
	cerr<<"time\t"<<(time(NULL)-t0)/60.0<<" min\n";
	return	0;
}

