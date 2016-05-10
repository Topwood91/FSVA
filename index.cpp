#include	<algorithm>
#include	<iostream>
#include	<stdint.h>
#include	<cstdlib>
#include	<cstring>
#include	<fstream>
#include	<sstream>
#include	<vector>
#include	<cmath>
#include	<map>
using	namespace	std;
const	unsigned	mod=(1UL<<29)-1;
const	unsigned	step=8;
unsigned	kmer;

struct	Data{
	uint32_t	key,	pos;
	bool	operator()(Data	X,	Data	Y){	return	X.key==Y.key?X.pos<Y.pos:X.key<Y.key;	}
};

class	Index{
private:
	string	ref;
public:
	bool	load_ref(const	char	*F);
	bool	make_index(const	char	*F);
};

bool	Index::load_ref(const	char	*F){
	char	code[256],	buf[65536];
	for(size_t	i=0;	i<256;	i++)	code[i]=4;
	code['A']=code['a']=0;	code['C']=code['c']=1;	code['G']=code['g']=2;	code['T']=code['t']=3;
	FILE	*f=fopen(F,	"rb");
	if(f==NULL)	return	false;
	fseek(f,	0,	SEEK_END);
	ref.reserve(ftell(f)+1);
	fclose(f);
	f=fopen(F,	"rt");
	if(f==NULL)	return	false;
	while(fgets(buf,	65536,	f)!=NULL){
		if(buf[0]=='>')	continue;
		for(char	*p=buf;	*p;	p++)	if(*p>=33)	ref.push_back(*(code+*p));
	}
	fclose(f);
	cerr<<"genome\t"<<ref.size()<<'\n';
	return	true;
}

bool	Index::make_index(const	char	*F){
	cerr<<"hashing\n";
	vector<Data>	data;	data.reserve(ref.size()/step+1);
	Data	d;
	for(size_t	i=0;	i<=ref.size()-kmer;	i+=step){
		uint64_t	h=0;	bool	hasn=false;
		for(size_t	j=0;	j<kmer;	j++){
			if(ref[i+j]==4){	hasn=true;	break;	}
			else	h=(h<<2)+ref[i+j];
		}
		if(hasn)	continue;
		d.key=h%mod;	d.pos=i;	data.push_back(d);
	}
	cerr<<"hash\t"<<data.size()<<'\n';
	cerr<<"sorting\n";
	sort(data.begin(),	data.end(),	Data());

	map<unsigned,	unsigned>	m;
	cerr<<"writing\n";
	string	fn=F;	fn+=".hash";
	ofstream	fo(fn.c_str(),	ios::binary);
	uint32_t	offset=data.size();
	fo.write((char*)&offset,	4);
	for(size_t	i=0;	i<data.size();	i++)	fo.write((char*)&data[i].pos,	4);
	size_t	last_key=0;
	for(size_t	i=0;	i<data.size();	){
		size_t	h=data[i].key,	n;	offset=i;
		for(size_t	j=last_key;	j<=h;	j++){	fo.write((char*)&offset,	4);	}
		last_key=h+1;
		for(n=i+1;	n<data.size()&&data[n].key==h;	n++);
		m[n-i]++;
		i=n;
	}
	offset=data.size();
	for(size_t	j=last_key;	j<=mod;	j++)	fo.write((char*)&offset,	4);
	fo.close();

	//fo.open("distribution.txt");
	//for(map<unsigned,	unsigned>::iterator	mi=m.begin();	mi!=m.end();	++mi)	fo<<mi->first<<'\t'<<mi->second<<'\n';
	//fo.close();
	return	true;
}

int	main(int	ac,	char	**av){
    kmer=31;
	if(ac<2){	cerr<<"index [options] <ref.fa>\n";
        cerr<<"options:\n";
        cerr<<"\t-r INT length of read. Using this argument, length of seed\n\t       can be calculated by program automatically. [150]\n";
        cerr<<"\t-l INT length of seed. [31]\n\n";
        cerr<<"\tDon't use both -r and -l\n\n";
        return	0;
    }
    unsigned readLen=0,kmer_temp=0;
    for(int it=1;it<ac;it++)
    {
        if(strcmp(av[it],"-l")==0)kmer_temp=atoi(av[it+1]);
        else if(strcmp(av[it],"-r")==0)readLen=atoi(av[it+1]);
    }
    if(readLen>=175)kmer=32;
    else if(readLen>=150)kmer=31;
    else if(readLen>=125)kmer=30;
    else if(readLen>=75)kmer=21;
    else if(readLen>0)kmer=16;
    if(kmer_temp!=0)kmer=kmer_temp;
	Index	i;
	if(!i.load_ref(av[ac-1]))	return	0;
	if(!i.make_index(av[ac-1]))	return	0;
	return	0;
}
