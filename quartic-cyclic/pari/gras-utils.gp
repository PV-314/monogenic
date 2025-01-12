isMOK(m)={
	return(1);
	if (m%9==0,return(0));
	if (m%25==0,return(0));
	if (m%49==0,return(0));
	if (m%2==0,
		if (m%16!=8,
			return(0);
		);
	);
	return(1);
}

areGAndMOK(g, m)={
	my(g2,gTmp,isProp1OK);
	
	return(1);
	if (g%9==0,return(0));
	if (g%25==0,return(0));
	if (g%49==0,return(0));

	\\ next we apply Gras' Proposition 1 in her PMB 1981 paper
	g2=g*g;
	if ((g2+4)%m!=0 && (g2-4)%m!=0,
		return(0);
	);
	\\ need more complicated logic here for case when m=8,
	\\ when both g^2+4 and g^2-4 can be divisible by 8
	isProp1OK=0;
	if ((g2+4)%m==0,
		gTmp=(g2+4)/m;
		if (issquare(gTmp),
			isProp1OK=1;
		);
	);
	if (isProp1OK==0,
		if ((g2-4)%m==0,
			gTmp=(g2-4)/m;
			if (issquare(gTmp),
				isProp1OK=1;
			);
		);
	);
	if (isProp1OK==0,
		return(0);
	);

	if (m%2==1,
		if(g%2==1,
			return(1);
		);
		if (g%8==4,
			return(1);
		);
		if (g%16==8,
			return(1);
		);
		return(0);
	);
	if (m%2==0,
		if (g%4==2 && m%16==8,
			return(1);
		);
		return(10);
	);
}

fullAreGAndMOK(g, m,dbg=0)={
	my(g2,gTmp,isProp1OK);
	
	if (g%9==0,return(0));
	if (g%25==0,return(0));
	if (g%49==0,return(0));
	if (m%9==0,return(0));
	if (m%25==0,return(0));
	if (m%49==0,return(0));

	\\ next we apply Gras' Proposition 1 in her PMB 1981 paper
	g2=g*g;
	if ((g2+4)%m!=0 && (g2-4)%m!=0,
		if(dbg!=0,
			print("   fullAreGAndMOK() failing for g=",g,", m=",m,", as (g2+4)%m=",(g2+4)%m," and (g2-4)%m=",(g2-4)%m);
		);
		return(0);
	);
	\\ need more complicated logic here for case when m=8,
	\\ when both g^2+4 and g^2-4 can be divisible by 8
	isProp1OK=0;
	if ((g2+4)%m==0,
		gTmp=(g2+4)/m;
		if (issquare(gTmp),
			isProp1OK=1;
		);
	);
	if (isProp1OK==0,
		if ((g2-4)%m==0,
			gTmp=(g2-4)/m;
			if (issquare(gTmp),
				isProp1OK=1;
			);
		);
	);
	if (isProp1OK==0,
		if(dbg!=0,
			print("   fullAreGAndMOK() failing Gras' Prop 1 for g=",g,", m=",m);
		);
		return(0);
	);

	if (m%2==1,
		if(g%2==1,
			return(1);
		);
		if (g%8==4,
			return(1);
		);
		if (g%16==8,
			return(1);
		);
		if(dbg!=0,
			print("   fullAreGAndMOK() failing odd m tests for g=",g,"=",g%16," mod 16, m=",m);
		);
	);
	if (m%2==0,
		if (g%4==2 && m%16==8,
			return(1);
		);
		if(dbg!=0,
			print("   fullAreGAndMOK() failing even m tests for g=",g,"=",g%4," mod 4, m=",m,"=",m%16," mod 16");
		);
	);
}

\\ use t for-loop rather than old calc_t() function due to concerns about how to determine t from Gras paper
\\ 17 Dec 2024
check_poly(a,b,g,x,y,z,dbg=0)={
	my(f,isFound,m,nf);
	
	m=a*a+b*b;
	printf("a=%5d, b=%6d, g=%10d, m=%10d, x=%5d, y=%5d, z=%5d, fullAreGAndMOK(g,m)=%1d, factor(g)=%s, factor(m)=%s\n",a,b,g,m,x,y,z,fullAreGAndMOK(g,m),factor(g),factor(m));
	isFound=0;
	for(t=0,10,
		if(isFound==0,
			f=calc_poly(a,b,g,t,x,y,z,dbg);
			print("   t=",t,", f=",f);
			if(polisirreducible(f) && denominator(content(f))==1,
				nf=bnfinit(f,1);
				print("   poldisc(f)=",poldisc(f),", nf.disc=",nf.disc);
				if(poldisc(f)==nf.disc,
					print("   MONOGENIC, poldisc(f)=",poldisc(f),", nf.disc=",nf.disc);
					isFound=1;
				);
				if(poldisc(f)!=nf.disc,
					print("   not monogenic, poldisc(f)/nf.disc=",poldisc(f)/nf.disc,", m=",factor(m));
				);
			);
			\\if(!polisirreducible(f),
			\\	print("REDUCIBLE: for z=",z,", f=",factor(f));
			\\);
		);
	);
}

\\ do a small search for t and print the polynomials found (think such values of t follow a mod 8 pattern)
\\ 18 Dec 2024
simple_check_poly(a,b,g,x,y,z,dbg=0)={
	my(d,f,isFound,m,nf);
	
	m=a*a+b*b;
	d=m*m*m*g*g; \\ discriminant of field (if m and g are "good")
	printf("   a=%5d, b=%6d, g=%10d, m=%10d, x=%5d, y=%5d, z=%5d, fullAreGAndMOK(g,m)=%1d, factor(g)=%s, factor(m)=%s\n",a,b,g,m,x,y,z,fullAreGAndMOK(g,m),factor(g),factor(m));
	isFound=0;
	for(t=0,16,
		f=calc_poly(a,b,g,t,x,y,z,dbg);
		if(polisirreducible(f) && denominator(content(f))==1,
			print("      t=",t,", f=",f);
			if(poldisc(f)==d,
				\\if(isFound==0,
				\\	nf=bnfinit(f,1);
				\\	if(poldisc(f)==nf.disc,
				\\		print("      MONOGENIC, poldisc(f)=",poldisc(f),", nf.disc=",nf.disc);
				\\	);
				\\	if(poldisc(f)!=nf.disc,
				\\		print("      not monogenic, poldisc(f)/nf.disc=",poldisc(f)/nf.disc,", m=",factor(m));
				\\	);
				\\);
			);
			isFound=1;
		);
	);
	return(isFound);
}

\\ this 
\\ 17 Dec 2024
full_check_poly(a,b,g,x,y,z,dbg=0)={
	my(f,m,nf,t);
	
	m=a*a+b*b;
	print("\n\na=",a,", b=",b,", g=",g,"=",factor(g),", m=",m,"=",factor(m),", x=",x,", y=",y,", z=",z);
	t=calc_t(g,m,x,y,z);
	cond=m*g;
	tau2=g*(a+b*I)*sqrt(m);
	tauBar2=g*(a-b*I)*sqrt(m);
	tau=sqrt(tau2);
	al=t+z*sqrt(m)+(x+y*I)*tau+(x-y*I)*conj(tau);
	al=real(al)/4;
	print("al=",al);
	f=calc_poly(a,b,g,t,x,y,z,dbg);
	print("f(al)=",subst(f,X,al));
	myPsi1=sqrt(g*sqrt(m)*(sqrt(m)+a)/2);
	myPsi2=-sqrt(g*sqrt(m)*(sqrt(m)+a)/2);
	myPsi3=sqrt(-g*sqrt(m)*(-sqrt(m)+a)/2);
	myPsi4=-sqrt(-g*sqrt(m)*(-sqrt(m)+a)/2);
	
	\\ using formula on page 3, line -2 of Gras' 1978 paper
	sigmaPsi=-b*g*sqrt(m)/2/myPsi1;
	th1=(t+z*sqrt(m)+2*x*myPsi1+2*y*myPsi2)/4;
	th2=(t+z*sqrt(m)+2*x*myPsi1+2*y*myPsi3)/4;
	th3=(t+z*sqrt(m)+2*x*myPsi1+2*y*myPsi4)/4;
	print("th1=",th1);
	print("th2=",th2);
	print("th3=",th3);
	thFromSigma=(t+z*sqrt(m)+2*x*myPsi1+2*y*sigmaPsi)/4;
	print("f(th1)=",subst(f,X,th1));
	print("f(th2)=",subst(f,X,th2));
	print("f(th3)=",subst(f,X,th3));
	print("f(thFromSigma)=",subst(f,X,thFromSigma));
	
	print("   t=",t,", f=",f);
	myDisc=m*m*m*g*g;
	print("poldisc(f)=",poldisc(f),", myDisc=",myDisc);
	if(denominator(content(f))==1,
		nf=bnfinit(f,1);
		print("   poldisc(f)=",poldisc(f),", nf.disc=",nf.disc);
		if(poldisc(f)==nf.disc,
			print("   MONOGENIC, poldisc(f)=",poldisc(f),", nf.disc=",nf.disc);
		);
		if(poldisc(f)!=nf.disc,
			print("   not monogenic, poldisc(f)/nf.disc=",poldisc(f)/nf.disc,", m=",factor(m));
		);
	);
}

\\ 15 Dec 2024
calc_poly(av,bv,gv,tv,xv,yv,zv,dbg=0)={
	my(bigS1,bigS2,bigS3,bigS4,f,fv,mv);
	
	if(dbg!=0,
		print("in calc_poly(): a=",av,", b=",bv,", g=",gv,", t=",tv,", x=",xv);
	);
	mv=av*av+bv*bv;
	fv=mv*gv;
	bigS1=tv;
	bigS2=tv*tv+mv*zv*zv-2*(xv*xv+yv*yv)*fv+2*(tv*tv-mv*zv*zv);
	bigS2=bigS2/8;
	
	bigS3=(tv*tv+mv*zv*zv-2*(xv*xv+yv*yv)*fv)*tv-2*mv*zv*(tv*zv-gv*(av*(xv*xv-yv*yv)-2*bv*xv*yv));
	bigS3=bigS3/16;

	bigS4=(tv*tv+mv*zv*zv-2*(xv*xv+yv*yv)*fv)*(tv*tv+mv*zv*zv-2*(xv*xv+yv*yv)*fv)-4*mv*(tv*zv-gv*(av*(xv*xv-yv*yv)-2*bv*xv*yv))*(tv*zv-gv*(av*(xv*xv-yv*yv)-2*bv*xv*yv));
	bigS4=bigS4/256;
	if(dbg!=0,
		print("in calc_poly(): bigS2=",bigS2);
		print("in calc_poly(): bigS3=",bigS3);
		print("in calc_poly(): bigS4=",bigS4);
	);
	
	f=X^4-bigS1*X^3+bigS2*X^2-bigS3*X+bigS4;
	return(f);
}

\\ uses equation (0) in Gras 1981 paper (actually from Hasse), which is not right
\\ 17 Dec 2024
DONOTUSE_calc_t(g,m,x,y,z)={
	my(t,v1,v2);
	
	if(m%2==1,
		t=z;
		v1=(t+z)/2;
		v2=(t-z)/2;
		if(v1%2!=(g*x)%2,
			t=z+2;
			v1=(t+z)/2;
			v2=(t-z)/2;
			if(v1%2!=(g*x)%2,
				print("BAD t for m odd: fails gx condition, g=",g,", m=",m,", t=",t,", x=",x,", y=",y,", z=",z,", (t+z)/2=",v1,", g*x=",g*x);
				return();
			);
		);
		if(v2%2!=(g*y)%2,
			print("BAD t for m odd: fails gy condition, g=",g,", m=",m,", t=",t,", x=",x,", y=",y,", z=",z,", (t-z)/2=",v2,", g*y=",g*y);
			return();
		);
	);
	if(m%2==0,
		if(z%2!=0,
			print("BAD z=",z,", for m even");
			return();
		);
		t=0;
	);
	return(t);
}