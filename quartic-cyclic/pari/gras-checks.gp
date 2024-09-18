\\ Gras Besancon 1981 paper Theoreme 1
\\ main entry point to other functions
\\ check_family1() at the end may also be of interest
\\ 17 May 2024
thm1_check(dbg=0)={
	my(bUB,condSet,m,v1,v2,y1,y2);

	condSet=Set();
	for(a=2,500,
		bUB=22200;
		if(a==1,bUB=3200);
		forstep(b=2,bUB,2,
			m=a*a+b*b;
			for(x=1,200,
				\\ from f1=b*(x*x-y*y)+2*a*x*y-(+2);
				v1=a*a*x*x+b*b*x*x-2*b;
				if(issquare(v1),
					v1s=sqrtint(v1);
					y1=(2*a*x+2*v1s)/2/b;
					if(denominator(y1)==1,
						condSet=check_eqn2(a,b,m,x,y1,condSet);
					);
					y2=(2*a*x-2*v1s)/2/b;
					if(denominator(y2)==1,
						condSet=check_eqn2(a,b,m,x,y2,condSet);
					);
				);

				\\ from f2=b*(x*x-y*y)+2*a*x*y-(-2);
				v2=a*a*x*x+b*b*x*x+2*b;
				if(issquare(v2),
					v2s=sqrtint(v2);
					y1=(2*a*x+2*v2s)/2/b;
					if(denominator(y1)==1,
						condSet=check_eqn2(a,b,m,x,y1,condSet);
					);
					y2=(2*a*x-2*v2s)/2/b;
					if(denominator(y2)==1,
						condSet=check_eqn2(a,b,m,x,y2,condSet);
					);
				);
			);
		);
		\\print("\n",condSet);
	);
	check_conductors(condSet);
}

\\ search for solutions when x=1 and y=2
\\ no further examples besides a=1, b=2,... and
\\ FOUND        : for a=   5, b=   6, g= 1523, m=      61, x=   1, y=   2, z=  85
\\ found
\\ 20 May 2024
thm1_check2(dbg=0)={
	my(b,b1,m,x,y);

	x=1;
	y=2;
	for(a=1,5000,
		if(a%20==0,
			print("starting a=",a);
		);
		b1=4*a-2;
		if(b1%6==0,
			b=b1/3;
			m=a*a+b*b;
			check_eqn2(a,b,m,x,y,condSet);
		);
		b1=4*a+2;
		if(b1%6==0,
			b=b1/3;
			m=a*a+b*b;
			check_eqn2(a,b,m,x,y,condSet);
		);
	);
}

\\ looking for real fields, so chi(-1)=+1
\\ eqn2 is the second equation in equation (1) in Gras' paper (part of her Theoreme 1)
check_eqn2(a,b,m,x,y,condSet)={
	my(gUB,isOK,v1,v2,v2Sqr,v3,v4,v4Sqr,v5);
	
	\\print("checking a=",a,", b=",b,", m=",m,", x=",x,", y=",y);
	gUB=10*1000*1000;
	for(g=1,gUB,
		isOK=1;
		if(m%2==1,
			if(g%4==2 || g%16==0,
				isOK=0;
			);
		);
		if(m%2==0,
			if(m%16!=8,
				isOK=0;
			);
			if(g%4!=2,
				isOK=0;
			);
		);
		if(isOK==1,
			\\ from m*(z^2-g*(x^2+y^2))^2-4*g*g-(+16)
			if(issquare(m*g*g+4*m),
				v1=sqrtint(m*g*g+4*m);
				v2Sqr=m*(m*g*x*x+m*g*y*y+2*v1);
				if(issquare(v2Sqr),
					v2=sqrtint(v2Sqr)/m;
					v3=-v2;
					if(denominator(v2)==1,
						printf("SOLUTION     : a=%4d, b=%4d, m=%6d, g=%7d, x=%4d, y=%4d, z=%6d, %6d\n",a,b,m,g,x,y,v2,v3);
						check_soln(a,b,g,m,x,y,v2);
						check_soln(a,b,g,m,x,y,v3);
						\\condSet=setunion(condSet,Set([m*g]));
					);
				);
				v4Sqr=m*(m*g*x*x+m*g*y*y-2*v1);
				if(issquare(v4Sqr),
					v4=sqrtint(v4Sqr)/m;
					v5=-v4;
					if(denominator(v4)==1,
						printf("\nSOLUTION     : a=%4d, b=%4d, m=%6d, g=%7d, x=%4d, y=%4d, z=%6d, %6d\n",a,b,m,g,x,y,v4,v5);
						check_soln(a,b,g,m,x,y,v4);
						check_soln(a,b,g,m,x,y,v5);
						\\condSet=setunion(condSet,Set([m*g]));
					);
				);
			);
			\\ from m*(z^2-g*(x^2+y^2))^2-4*g*g-(-16)
			if(issquare(m*g*g-4*m),
				v1=sqrtint(m*g*g-4*m);
				v2Sqr=m*(m*g*x*x+m*g*y*y+2*v1);
				if(issquare(v2Sqr),
					v2=sqrtint(v2Sqr)/m;
					v3=-v2;
					if(denominator(v2)==1,
						printf("\nSOLUTION     : a=%4d, b=%4d, m=%6d, g=%7d, x=%4d, y=%4d, z=%6d, %6d\n",a,b,m,g,x,y,v2,v3);
						check_soln(a,b,g,m,x,y,v2);
						check_soln(a,b,g,m,x,y,v3);
						\\condSet=setunion(condSet,Set([m*g]));
					);
				);
				v4Sqr=m*(m*g*x*x+m*g*y*y-2*v1);
				if(issquare(v4Sqr),
					v4=sqrtint(v4Sqr)/m;
					v5=-v4;
					if(denominator(v4)==1,
						printf("\nSOLUTION     : a=%4d, b=%4d, m=%6d, g=%7d, x=%4d, y=%4d, z=%6d, %6d\n",a,b,m,g,x,y,v4,v5);
						check_soln(a,b,g,m,x,y,v4);
						check_soln(a,b,g,m,x,y,v5);
						\\condSet=setunion(condSet,Set([m*g]));
					);
				);
			);
		);
	);
	return(condSet);
}

\\ try to get a monogenic field from a solution obtained to the two equations
\\ in equation (1) within Theoreme 1 of Gras' 1981 paper
\\ 17 May 2024
check_soln(a,b,g,m,x,y,z,dbg=0)={
	my(condList,f,f1,f2,g1,g2,isCondOK,myPsi1,myPsi2,myPsi3,myPsi4,nf,nf1,nf2,psiList,th1,th2,th3,th4,thPolySet,tv);
	
	\\ turn off the conductor check as it is too expensive time-wise as m*g grows
	useCondCheck=0;
	\\if(!issquarefree(m*g/gcd(m*g,65536*65536)),
	\\	printf("NOT SQRFREE  : a=%4d, b=%4d, m=%6d, g=%7d, x=%4d, y=%4d, z=%6d\n",a,sqrtint(m-a*a),m,g,x,y,z);
	\\	return();
	\\);

	if(useCondCheck!=0,
		if(dbg!=0,	
			print("starting check_soln() with a=",a,", b=",b,", g=",g,", m=",m);
		);
		condList=get_conductor(m*g,dbg);
		if(dbg!=0,	
			print("done get_conductor() in check_soln() with a=",a,", b=",b,", g=",g,", m=",m);
		);
	);
	
	myPsi1=sqrt(g*sqrt(m)*(sqrt(m)+a)/2);
	myPsi2=-sqrt(g*sqrt(m)*(sqrt(m)+a)/2);
	myPsi3=sqrt(-g*sqrt(m)*(-sqrt(m)+a)/2);
	myPsi4=-sqrt(-g*sqrt(m)*(-sqrt(m)+a)/2);
	psiList=List([myPsi1,myPsi2,myPsi3,myPsi4]);
	f1=(x1-myPsi1)*(x1-myPsi2)*(x1-myPsi3)*(x1-myPsi4);
	g1=get_monic_poly(f1);
	nf1=bnfinit(g1,1);

	if(useCondCheck!=0,
		isCondOK=0;
		for(i=1,length(condList),
			f=condList[i];
			nf=bnfinit(f,1);
			if(nf1.disc==nf.disc,
				isCondOK=1;
			);
		);
		if(isCondOK==0,
			printf("NOT CONDUCTOR: a=%4d, b=%4d, m=%6d, g=%7d, x=%4d, y=%4d, z=%6d\n",a,sqrtint(m-a*a),m,g,x,y,z);
			return();
		);
		printf("CONDUCTOR    : a=%4d, b=%4d, m=%6d, g=%7d, x=%4d, y=%4d, z=%6d\n",a,sqrtint(m-a*a),m,g,x,y,z);
	);

	if(m%2==0 && z%2==0,
		tv=4;
		thPolySet=get_theta_polys(a,b,g,m,tv,x,y,z,psiList);
		for(i=1,length(thPolySet),
			g2=thPolySet[i];
			if(!polisirreducible(g2),
				printf("REDUCIBLE    : a=%4d, b=%4d, g=%5d, m=%8d, x=%4d, y=%4d, z=%4d\n",a,b,g,m,x,y,z);
				print("tried tv=",tv);
				print("g2=",g2);
				1/0;
			);
			if(polisirreducible(g2),
				nf2=bnfinit(g2,1);
				if(dbg!=0 && poldisc(g2)!=nf2.disc,
					printf("for a=%4d, b=%4d, g=%5d, m=%8d, x=%4d, y=%4d, z=%4d\n",a,b,g,m,x,y,z);
					print("tried tv=",tv);
					print("g1=",g1);
					printf("BAD DISC: disc(g2)=%20d, nf2.disc=%20d=%s, g2=%s, Gal=%s\n",poldisc(g2),nf2.disc,factor(nf2.disc),g2,polgalois(g2));
				);
				if(poldisc(g2)==nf2.disc,
					printf("\nFOUND        : a=%4d, b=%4d, g=%5d, m=%8d, x=%4d, y=%4d, z=%4d\n",a,b,g,m,x,y,z);
					print("using tv=",tv);
					print("g1=",g1);
					printf("disc(g2)=nf2.disc=%30d=%s, g2=%s, Gal=%s\n",poldisc(g2),factor(nf2.disc),g2,polgalois(g2));
				);
			);
		);
	);
	
	if(m%2==1,
		tv=z%2;
		if(dbg!=0,
			print("for m odd, starting with tv=",tv," and working forward with values of tv");
		);
		while(((tv+z)/2)%2!=(g*x)%2 || ((tv-z)/2)%2!=(g*y)%2,
			tv=tv+2;
		);
		if(dbg!=0,
			print("for m odd, found tv=",tv);
		);
		thPolySet=get_theta_polys(a,b,g,m,tv,x,y,z,psiList);
		for(i=1,length(thPolySet),
			g2=thPolySet[i];
			if(!polisirreducible(g2),
				printf("REDUCIBLE    : a=%4d, b=%4d, g=%5d, m=%8d, x=%4d, y=%4d, z=%4d\n",a,b,g,m,x,y,z);
				print("tried tv=",tv);
				print("g2=",g2);
			);
			if(polisirreducible(g2),
				nf2=bnfinit(g2,1);
				if(dbg!=0 && poldisc(g2)!=nf2.disc,
					printf("for a=%4d, b=%4d, g=%5d, m=%8d, x=%4d, y=%4d, z=%4d\n",a,b,g,m,x,y,z);
					print("tried tv=",tv);
					print("g1=",g1);
					printf("BAD DISC: disc(g2)=%20d, nf2.disc=%20d=%s, g2=%s, Gal=%s\n",poldisc(g2),nf2.disc,factor(nf2.disc),g2,polgalois(g2));
				);
				if(poldisc(g2)==nf2.disc,
					printf("\nFOUND        : a=%4d, b=%4d, g=%5d, m=%8d, x=%4d, y=%4d, z=%4d\n",a,b,g,m,x,y,z);
					print("using tv=",tv);
					print("g1=",g1);
					printf("disc(g2)=nf2.disc=%30d=%s, g2=%s, Gal=%s\n",poldisc(g2),factor(nf2.disc),g2,polgalois(g2));
				);
			);
		);

		tv=z%2-2;
		if(dbg!=0,
			print("for m odd, starting with tv=",tv," and working backward with values of tv");
		);
		while(((tv+z)/2)%2!=(g*x)%2 || ((tv-z)/2)%2!=(g*y)%2,
			tv=tv-2;
		);
		if(dbg!=0,
			print("for m odd, found tv=",tv);
		);
		thPolySet=get_theta_polys(a,b,g,m,tv,x,y,z,psiList);
		for(i=1,length(thPolySet),
			g2=thPolySet[i];
			if(!polisirreducible(g2),
				printf("REDUCIBLE    : a=%4d, b=%4d, g=%5d, m=%8d, x=%4d, y=%4d, z=%4d\n",a,b,g,m,x,y,z);
				print("tried tv=",tv);
				print("g2=",g2);
			);
			if(polisirreducible(g2),
				nf2=bnfinit(g2,1);
				if(dbg!=0 && poldisc(g2)!=nf2.disc,
					printf("for a=%4d, b=%4d, g=%5d, m=%8d, x=%4d, y=%4d, z=%4d\n",a,b,g,m,x,y,z);
					print("tried tv=",tv);
					print("g1=",g1);
					printf("BAD DISC: disc(g2)=%20d, nf2.disc=%20d=%s, g2=%s, Gal=%s\n",poldisc(g2),nf2.disc,factor(nf2.disc),g2,polgalois(g2));
				);
				if(poldisc(g2)==nf2.disc,
					printf("\nFOUND        : a=%4d, b=%4d, g=%5d, m=%8d, x=%4d, y=%4d, z=%4d\n",a,b,g,m,x,y,z);
					print("using tv=",tv);
					print("g1=",g1);
					printf("disc(g2)=nf2.disc=%30d=%s, g2=%s, Gal=%s\n",poldisc(g2),factor(nf2.disc),g2,polgalois(g2));
				);
			);
		);
	);
}

\\ returns a set of all possible theta polynomials
\\ for a given tuple (a,b,g,m,t,x,y,z)
\\ 19 May 2024
get_theta_polys(av,bv,gv,mv,tv,xv,yv,zv,psiList)={
	my(fTh,pSet);
	
	pSet=Set();
	for(i1=1,4,
	for(i2=1,4,
	for(i3=1,4,
	for(i4=1,4,
		if(i1!=1 && i2!=2 && i3!=3 && i4!=4 && i1!=i2 && i1!=i3 && i1!=i4 && i2!=i3 && i2!=i4 && i3!=i4,
			fTh=get_single_theta_poly(av,bv,gv,mv,tv,xv,yv,zv,+1,+1,+1,+1,i1,i2,i3,i4,psiList);
			if(fTh!=0,
				pSet=setunion(pSet,Set([fTh]));
			);
			fTh=get_single_theta_poly(av,bv,gv,mv,tv,xv,yv,zv,+1,+1,+1,-1,i1,i2,i3,i4,psiList);
			if(fTh!=0,
				pSet=setunion(pSet,Set([fTh]));
			);
			fTh=get_single_theta_poly(av,bv,gv,mv,tv,xv,yv,zv,+1,+1,-1,+1,i1,i2,i3,i4,psiList);
			if(fTh!=0,
				pSet=setunion(pSet,Set([fTh]));
			);
			fTh=get_single_theta_poly(av,bv,gv,mv,tv,xv,yv,zv,+1,+1,-1,-1,i1,i2,i3,i4,psiList);
			if(fTh!=0,
				pSet=setunion(pSet,Set([fTh]));
			);
			fTh=get_single_theta_poly(av,bv,gv,mv,tv,xv,yv,zv,+1,-1,+1,+1,i1,i2,i3,i4,psiList);
			if(fTh!=0,
				pSet=setunion(pSet,Set([fTh]));
			);
			fTh=get_single_theta_poly(av,bv,gv,mv,tv,xv,yv,zv,+1,-1,+1,-1,i1,i2,i3,i4,psiList);
			if(fTh!=0,
				pSet=setunion(pSet,Set([fTh]));
			);
			fTh=get_single_theta_poly(av,bv,gv,mv,tv,xv,yv,zv,+1,-1,-1,+1,i1,i2,i3,i4,psiList);
			if(fTh!=0,
				pSet=setunion(pSet,Set([fTh]));
			);
			fTh=get_single_theta_poly(av,bv,gv,mv,tv,xv,yv,zv,+1,-1,-1,-1,i1,i2,i3,i4,psiList);
			if(fTh!=0,
				pSet=setunion(pSet,Set([fTh]));
			);
			
			fTh=get_single_theta_poly(av,bv,gv,mv,tv,xv,yv,zv,-1,+1,+1,+1,i1,i2,i3,i4,psiList);
			if(fTh!=0,
				pSet=setunion(pSet,Set([fTh]));
			);
			fTh=get_single_theta_poly(av,bv,gv,mv,tv,xv,yv,zv,-1,+1,+1,-1,i1,i2,i3,i4,psiList);
			if(fTh!=0,
				pSet=setunion(pSet,Set([fTh]));
			);
			fTh=get_single_theta_poly(av,bv,gv,mv,tv,xv,yv,zv,-1,+1,-1,+1,i1,i2,i3,i4,psiList);
			if(fTh!=0,
				pSet=setunion(pSet,Set([fTh]));
			);
			fTh=get_single_theta_poly(av,bv,gv,mv,tv,xv,yv,zv,-1,+1,-1,-1,i1,i2,i3,i4,psiList);
			if(fTh!=0,
				pSet=setunion(pSet,Set([fTh]));
			);
			fTh=get_single_theta_poly(av,bv,gv,mv,tv,xv,yv,zv,-1,-1,+1,+1,i1,i2,i3,i4,psiList);
			if(fTh!=0,
				pSet=setunion(pSet,Set([fTh]));
			);
			fTh=get_single_theta_poly(av,bv,gv,mv,tv,xv,yv,zv,-1,-1,+1,-1,i1,i2,i3,i4,psiList);
			if(fTh!=0,
				pSet=setunion(pSet,Set([fTh]));
			);
			fTh=get_single_theta_poly(av,bv,gv,mv,tv,xv,yv,zv,-1,-1,-1,+1,i1,i2,i3,i4,psiList);
			if(fTh!=0,
				pSet=setunion(pSet,Set([fTh]));
			);
			fTh=get_single_theta_poly(av,bv,gv,mv,tv,xv,yv,zv,-1,-1,-1,-1,i1,i2,i3,i4,psiList);
			if(fTh!=0,
				pSet=setunion(pSet,Set([fTh]));
			);
		);
	);
	);
	);
	);
	return(pSet);
}

\\ av, bv and gv are only here for logging:
\\ returns 0 in case of reducible or coefficients not near an integer
get_single_theta_poly(av,bv,gv,mv,tv,xv,yv,zv,s1,s2,s3,s4,i1,i2,i3,i4,psiList)={
	my(fTh,th1,th2,th3,th4);
	
	th1=(tv+s1*zv*sqrt(mv)+2*xv*psiList[1]+2*yv*psiList[i1])/4;
	th2=(tv+s2*zv*sqrt(mv)+2*xv*psiList[2]+2*yv*psiList[i2])/4;
	th3=(tv+s3*zv*sqrt(mv)+2*xv*psiList[3]+2*yv*psiList[i3])/4;
	th4=(tv+s4*zv*sqrt(mv)+2*xv*psiList[4]+2*yv*psiList[i4])/4;
	fTh=(x-th1)*(x-th2)*(x-th3)*(x-th4);
	fTh=get_monic_poly(fTh);
	if(fTh!=0 && !polisirreducible(fTh),
		\\print(mv,tv,xv,yv,zv,s1,s2,s3,s4,i1,i2,i3,i4);
		printf("REDUCIBLE    : fTh=%s, for a=%5d, b=%5d, g=%5d, m=%5d, t=%3d, x=%4d, y=%4d, z=%4d, s1=%2d, s2=%2d, s3=%2d, s4=%2d, i1=%2d, i2=%2d, i3=%2d, i4=%2d, psiList=%s\n",fTh,av,bv,gv,mv,tv,xv,yv,zv,s1,s2,s3,s4,i1,i2,i3,i4,psiList);
		return(0);
	);
	return(fTh);
}

\\ helper function
\\ has an error or returns 0 if there are problems
\\ only logs integer error in debug, as often (e.g., for use with get_single_theta_poly()),
\\ the coefficients may not be integers as it is not the right poly
\\ (recall that we do not know which one is right, so we test all possibilities)
\\ 19 May 2024
get_monic_poly(f2,dbg=0)={
	my(cf,g2,tol);
	
	tol=0.000001;
	g2=0;
	for(i=0,poldegree(f2),
		cf=polcoef(f2,i);
		if(imag(cf)>tol,
			print("BAD: coeff not real for x^",i," with f2=",f2);
			1/0;
		);
		cf=real(cf);
		if(abs(cf-round(cf))>tol,
			if(dbg!=0,
				print("coeff not near an integer for x^",i," with f2=",f2);
			);
			return(0);
		);
		cf=round(cf);
		g2=g2+cf*x^i;
	);
	return(g2);
}

\\ pulled out for testing
\\ 18 May 2024
get_conductor(n,dbg=0)={
	my(currSet,f,gal,p1Set,pSet);
	
	pSet=Set(List(polsubcyclo(n,4)));
	
	\\ now make sure that we get the actual conductor, by removing non-mins
	nDivs=divisors(n);
	for(j=1,length(nDivs)-1,
		currSet=Set(List(polsubcyclo(nDivs[j],4)));
		pSet=setminus(pSet,currSet);
	);
	pSet1=Set(); \\ filter out non-cyclic ones
	for(j=1,length(pSet),
		f=pSet[j];
		gal=polgalois(f);
		\\ check that it is cyclic (polsubcyclo(n,4) just brings back degree 4 fields, nothing says they have to be cyclic)
		if(gal[1]==4 && gal[2]==-1,
			pSet1=setunion(pSet1,Set([f]));
		);
	);
	if(dbg!=0,
		print("in get_conductor(): for conductor=",n,", found ",length(pSet1)," possible polys");
	);
	return(pSet1);
}

\\ copied from Pari user guide
\\ 17 May 2024
mysubcyclo(n, d = -1)={
	my(bnr,L,IndexBound);

	IndexBound = if (d < 0, n, [d]);
	bnr = bnrinit(bnfinit(y), [n,[1]]);
	L = subgrouplist(bnr, IndexBound, 1);
	vector(#L,i, galoissubcyclo(bnr,L[i]));
}

\\ 17 May 2024
check_family1()={
	for(z=1,310,
		f=x^4-z*x^3+(-3/8*z^6+2*z^4-37/8*z^2+4)*x^2+(19/16*z^7-9/2*z^5+135/16*z^3-6*z-1/8*z^9)*x+2+9/4*z^6-107/128*z^8+1/4*z^2+5/32*z^10-3/256*z^12-719/256*z^4;
		if(denominator(content(f))==1 && issquarefree(z^2-2),
			nf=bnfinit(f,1);
			if(poldisc(f)==nf.disc,
				print("GOOD: z=",z,", nf.disc=",nf.disc,", poldisc(f)=",poldisc(f),", Gal(f)=",polgalois(f));
			);
			if(poldisc(f)!=nf.disc,
				printf("BAD:  z=%3d, ratio=%10d, nf.disc=%20d, poldisc(f)=%20d, Gal(f)=%s\n",z,(poldisc(f)/nf.disc),nf.disc,poldisc(f),polgalois(f));
			);
		);
	);
}
