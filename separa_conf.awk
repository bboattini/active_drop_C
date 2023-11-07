#!/usr/bin/awk -f

BEGIN{
    L=240;
    i=0;
    linha=0;
    count=0;
    saida = sprintf("%s_%d",file,count);
    saida1 = sprintf("%s_yHalf_%d",file,count);
    print saida
}
{
    if($0 !~ /^#/)  {
	if($0 !~ /./)  {
	    linha++;
	    print linha,count;
	    if(linha%2==0) {
		count++;
		ind=linha/2;
		saida = sprintf("%s_%d",file,ind);
		saida1 = sprintf("%s_yHalf_%d",file,ind);
		count++;
	    #print count
	    }
	}
	else {
	    #if((count==0)||(count==1)||(count==2)||(count==3)||(count==4)||(count==5)||(count==6)) {
	    #print count
	    num=int($1);
	    ss=int($2);
	    x = num%L;
	    y = int(num/L)%L;
	    z = int(int(num)/(int(L*L)));
	    s = x + y*L + z*L*L;
	    print $1,x,y,z,ss,s,$1-s >> saida
	    if(y==120) {
		print $1,x,y,z,ss,s,$1-s >> saida1
	    }
	    #}
	}
	i++;	
	#print i;
    }
}
END{
    print linha,i;
}
