void	radix_sort(unsigned	*begin,	unsigned *end){
	unsigned	*begin1=new	unsigned[end-begin],	*end1=begin1+(end-begin),	*temp;
    for(unsigned	shift=0;	shift<32;	shift+=8){
    	unsigned	count[0x100]={};
        for(unsigned	*p=begin;	p!=end;	p++)	count[(*p>>shift)&0xFF]++;
		unsigned	*bucket[0x100],	*q=begin1;
        for(unsigned	i=0;	i<0x100;	q+=count[i++])	bucket[i]=q;
        for(unsigned	*p=begin;	p!=end;	p++)	*bucket[(*p>>shift)&0xFF]++=*p;
		temp=begin;	begin=begin1;	begin1=temp;
		temp=end;	end=end1;	end1=temp;
    }
	delete[] begin1;
}
