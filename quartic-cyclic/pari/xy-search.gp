read("gras-utils.gp");

\\ for i=1,
\\ it seems that the constant coefficient of the f3 polynomial
\\ (j^2+1) * ( j^3+3*x-2 ) (here j is used as y)
\\ from interp([8,9,10,11,12,13,14],[34710,61828,103828,166164,255490,379780,548448],x);factor(%);
\\ the leading coefficient of the f3 polynomial is
\\ (x^2+y^2)^3
\\ 26 Dec 2024
general_search_all(dbg=0)={
	for(i=2,31,
	for(j=i+1,i+30,
		print("\nstarting x=",i,", y=",j);
		if(gcd(i,j)==1,
			general_xy_search(i,j);
		);
	);
	);
	for(i=1,31,
	for(j=i+31,i+50,
		print("\nstarting x=",i,", y=",j);
		if(gcd(i,j)==1,
			general_xy_search(i,j);
		);
	);
	);
}

\\ restart:b1:=expand(solve(b*(x^2-y^2)+2*a*x*y=2,b));b1a:=expand(subs(a=(y^2-x^2)*a1+z,b1));factor(coeff(%,a1,1));simplify(expand(b1a-(2*x*y*a1)));
\\ b2:=expand(solve(b*(x^2-y^2)+2*a*x*y=-2,b));b2a:=expand(subs(a=(y^2-x^2)*a1+z,b2)):factor(coeff(%,a1,1));simplify(expand(b2a-(2*x*y*a1)));
\\ 25 Dec 2024
general_xy_search(x,y,dbg=0)={
	my(a1,a10,a11,a2,a20,a21,b1,b10Coeff,b2,b20Coeff,c,eqn1V,eqn2aV,eqn2bV,g1a,g1b,g2a,g2b,m1,m2,v1,v2,
	id1a,id1b,id2a,id2b,id3a,id3b,id4a,id4b,id5a,id5b,id6a,id6b,id7a,id7b,id8a,id8b,
	x1aList,x1bList,x2aList,x2bList,x3aList,x3bList,x4aList,x4bList,x5aList,x5bList,x6aList,x6bList,x7aList,x7bList,x8aList,x8bList,
	y1aList,y1bList,y2aList,y2bList,y3aList,y3bList,y4aList,y4bList,y5aList,y5bList,y6aList,y6bList,y7aList,y7bList,y8aList,y8bList);
	
	if(y<=x,
		print("we assume that x<y");
		return();
	);
	if(gcd(x,y)>1,
		print("gcd(x,y)=",gcd(x,y),">1, so no solution to the first equation in Gras' criterion");
		return();
	);
	
	a11=y^2-x^2;
	if(a11%2==0,
		a11=a11/2;
	);
	a10=-1;
	a21=a11;
	a20=-1;
	for(z=0,a11-1,
		\\ eq'n=2 case: restart:b1:=expand(solve(b*(x^2-y^2)+2*a*x*y=2,b));b1a:=expand(subs(a=(y^2-x^2)*a1+z,b1));factor(coeff(%,a1,1));simplify(expand(b1a-(2*x*y*a1)));
		if(a10==-1 && (2*(x*y*z-1))%(y^2-x^2)==0,
			a10=z;
		);
		\\ eq'n=-2 case: b2:=expand(solve(b*(x^2-y^2)+2*a*x*y=-2,b));b2a:=expand(subs(a=(y^2-x^2)*a1+z,b2)):factor(coeff(%,a1,1));simplify(expand(b2a-(2*x*y*a1)));
		if(a20==-1 && (2*(x*y*z+1))%(y^2-x^2)==0,
			a20=z;
		);
	);
	print("a11=",a11,", a10=",a10,": a21=",a21,", a20=",a20);
	if(a10==-1 || a20==-1,
		print("BAD: could not solve for integer a's and b's");
		return();
	);
	b10Coeff=-2/(y^2-x^2);
	b20Coeff=2/(y^2-x^2);

	id1a="G1A-P2P16";
	id1b="G1A-M2P16";
	id2a="G1A-P2M16";
	id2b="G1A-M2M16";
	id3a="G1B-P2P16";
	id3b="G1B-M2P16";
	id4a="G1B-P2M16";
	id4b="G1B-M2M16";

	id5a="G2A-P2P16";
	id5b="G2A-M2P16";
	id6a="G2A-P2M16";
	id6b="G2A-M2M16";
	id7a="G2B-P2P16";
	id7b="G2B-M2P16";
	id8a="G2B-P2M16";
	id8b="G2B-M2M16";
	
	x1aList=List();
	y1aList=List();
	x1bList=List();
	y1bList=List();
	x2aList=List();
	y2aList=List();
	x2bList=List();
	y2bList=List();
	x3aList=List();
	y3aList=List();
	x3bList=List();
	y3bList=List();
	x4aList=List();
	y4aList=List();
	x4bList=List();
	y4bList=List();

	x5aList=List();
	y5aList=List();
	x5bList=List();
	y5bList=List();
	x6aList=List();
	y6aList=List();
	x6bList=List();
	y6bList=List();
	x7aList=List();
	y7aList=List();
	x7bList=List();
	y7bList=List();
	x8aList=List();
	y8aList=List();
	x8bList=List();
	y8bList=List();
		
	c=x*x+y*y;
	for(v=1,300000,
		if(dbg!=0 && v%10000==0,print("v=",v));
		\\ from 2 on RHS in first equation
		a1=a11*v+a10;
		b1=2*a1*x*y/(y^2-x^2)+b10Coeff;
		if(denominator(b1)!=1,
			print("BAD: b1 is not an integer for v=",v,", a1=",a1,", b1=",b1);
			return();
		);
		m1=a1*a1+b1*b1;
		
		\\ from -2 on RHS in first equation
		a2=a21*v+a20;
		b2=2*a2*x*y/(y^2-x^2)+2/(y^2-x^2);
		if(denominator(b2)!=1,
			print("BAD: b2 is not an integer for v=",v,", a2=",a2,", b2=",b2);
			return();
		);
		m2=a2*a2+b2*b2;

		v1=m1*c*c;
		if(issquare(v1+4),
			g1a=sqrtint(v1+4);
			eqn1V=b1*(x^2-y^2)+2*a1*x*y;
			if(abs(eqn1V)!=2,
				printf("BAD: (G1A), a=%7d, b=%7d, g=%10d, m=%14d, x=%2d, y=%2d, c=%6d, eqn1V=%2d, eqn2V=%s*%s\n",a1,b1,g1a,m1,x,y,c,eqn1V,content(eqn2aV),factor(eqn2aV));
				return();
			);
			eqn2aV=m1*(z*z-g1a*(x*x+y*y))*(z*z-g1a*(x*x+y*y))-4*g1a*g1a-16;
			if(!polisirreducible(eqn2aV),
				if(eqn1V==2,
					\\if(has_linear_factor(eqn2aV) || dbg!=0,
					if(dbg!=0,
						printf("(%s): a=%7d, b=%7d, g=%10d, m=%14d, x=%2d, y=%2d, c=%6d, eqn1V=%2d, eqn2V=%s*%s\n",id1a,a1,b1,g1a,m1,x,y,c,eqn1V,content(eqn2aV),factor(eqn2aV));
					);
					if(length(x1aList)<10,
						listput(x1aList,a1);
						listput(y1aList,c*g1a-2*c);
					);
				);
				if(eqn1V==-2,
					\\if(has_linear_factor(eqn2aV) || dbg!=0,
					if(dbg!=0,
						printf("(%s): a=%7d, b=%7d, g=%10d, m=%14d, x=%2d, y=%2d, c=%6d, eqn1V=%2d, eqn2V=%s*%s\n",id1b,a1,b1,g1a,m1,x,y,c,eqn1V,content(eqn2aV),factor(eqn2aV));
					);
					if(length(x1bList)<10,
						listput(x1bList,a1);
						listput(y1bList,c*g1a-2*c);
					);
				);
			);
			eqn2bV=m1*(z*z-g1a*(x*x+y*y))*(z*z-g1a*(x*x+y*y))-4*g1a*g1a+16;
			if(!polisirreducible(eqn2bV),
				if(eqn1V==2,
					\\if(has_linear_factor(eqn2bV) || dbg!=0,
					if(dbg!=0,
						printf("(%s): a=%7d, b=%7d, g=%10d, m=%14d, x=%2d, y=%2d, c=%6d, eqn1V=%2d, eqn2V=%s*%s\n",id2a,a1,b1,g1a,m1,x,y,c,eqn1V,content(eqn2bV),factor(eqn2bV));
					);
					if(length(x2aList)<10,
						listput(x2aList,a1);
						listput(y2aList,c*g1a-2*c);
					);
				);
				if(eqn1V==-2,
					\\if(has_linear_factor(eqn2bV) || dbg!=0,
					if(dbg!=0,
						printf("(%s): a=%7d, b=%7d, g=%10d, m=%14d, x=%2d, y=%2d, c=%6d, eqn1V=%2d, eqn2V=%s*%s\n",id2b,a1,b1,g1a,m1,x,y,c,eqn1V,content(eqn2bV),factor(eqn2bV));
					);
					if(length(x2bList)<10,
						listput(x2bList,a1);
						listput(y2bList,c*g1a-2*c);
					);
				);
			);
		);
		if(issquare(v1-4),
			g1b=sqrtint(v1-4);
			eqn1V=b1*(x^2-y^2)+2*a1*x*y;
			if(abs(eqn1V)!=2,
				printf("BAD: (G1B), a=%7d, b=%7d, g=%10d, m=%14d, x=%2d, y=%2d, c=%6d, eqn1V=%2d\n",a1,b1,g1b,m1,x,y,c,eqn1V);
				return();
			);
			eqn2aV=m1*(z*z-g1b*(x*x+y*y))*(z*z-g1b*(x*x+y*y))-4*g1b*g1b-16;
			if(!polisirreducible(eqn2aV),
				if(eqn1V==2,
					\\if(has_linear_factor(eqn2aV) || dbg!=0,
					if(dbg!=0,
						printf("(%s): a=%7d, b=%7d, g=%10d, m=%14d, x=%2d, y=%2d, c=%6d, eqn1V=%2d, eqn2V=%s*%s\n",id3a,a1,b1,g1b,m1,x,y,c,eqn1V,content(eqn2aV),factor(eqn2aV));
					);
					if(length(x3aList)<10,
						listput(x3aList,a1);
						listput(y3aList,c*g1b-2*c);
					);
				);
				if(eqn1V==-2,
					\\if(has_linear_factor(eqn2aV) || dbg!=0,
					if(dbg!=0,
						printf("(%s): a=%7d, b=%7d, g=%10d, m=%14d, x=%2d, y=%2d, c=%6d, eqn1V=%2d, eqn2V=%s*%s\n",id3b,a1,b1,g1b,m1,x,y,c,eqn1V,content(eqn2aV),factor(eqn2aV));
					);
					if(length(x3bList)<10,
						listput(x3bList,a1);
						listput(y3bList,c*g1b-2*c);
					);
				);
			);
			eqn2bV=m1*(z*z-g1b*(x*x+y*y))*(z*z-g1b*(x*x+y*y))-4*g1b*g1b+16;
			if(!polisirreducible(eqn2bV),
				if(eqn1V==2,
					\\if(has_linear_factor(eqn2bV) || dbg!=0,
					if(dbg!=0,
						printf("(%s): a=%7d, b=%7d, g=%10d, m=%14d, x=%2d, y=%2d, c=%6d, eqn1V=%2d, eqn2V=%s*%s\n",id4a,a1,b1,g1b,m1,x,y,c,eqn1V,content(eqn2bV),factor(eqn2bV));
					);
					if(length(x4aList)<10,
						listput(x4aList,a1);
						listput(y4aList,c*g1b-2*c);
					);
				);
				if(eqn1V==-2,
					\\if(has_linear_factor(eqn2bV) || dbg!=0,
					if(dbg!=0,
						printf("(%s): a=%7d, b=%7d, g=%10d, m=%14d, x=%2d, y=%2d, c=%6d, eqn1V=%2d, eqn2V=%s*%s\n",id4b,a1,b1,g1b,m1,x,y,c,eqn1V,content(eqn2bV),factor(eqn2bV));
					);
					if(length(x4bList)<10,
						listput(x4bList,a1);
						listput(y4bList,c*g1b-2*c);
					);
				);
			);
		);
		v2=m2*c*c;
		if(issquare(v2+4),
			g2a=sqrtint(v2+4);
			eqn1V=b2*(x^2-y^2)+2*a2*x*y;
			if(abs(eqn1V)!=2,
				printf("BAD: (G2A), a=%7d, b=%7d, g=%10d, m=%14d, x=%2d, y=%2d, c=%6d, eqn1V=%2d\n",a2,b2,g2a,m2,x,y,c,eqn1V);
				return();
			);
			eqn2aV=m2*(z*z-g2a*(x*x+y*y))*(z*z-g2a*(x*x+y*y))-4*g2a*g2a-16;
			if(!polisirreducible(eqn2aV),
				if(eqn1V==2,
					\\if(has_linear_factor(eqn2aV) || dbg!=0,
					if(dbg!=0,
						printf("(%s): a=%7d, b=%7d, g=%10d, m=%14d, x=%2d, y=%2d, c=%6d, eqn1V=%2d, eqn2V=%s*%s\n",id5a,a2,b2,g2a,m2,x,y,c,eqn1V,content(eqn2aV),factor(eqn2aV));
					);
					if(length(x5aList)<10,
						listput(x5aList,a2);
						listput(y5aList,c*g2a-2*c);
					);
				);
				if(eqn1V==-2,
					\\if(has_linear_factor(eqn2aV) || dbg!=0,
					if(dbg!=0,
						printf("(%s): a=%7d, b=%7d, g=%10d, m=%14d, x=%2d, y=%2d, c=%6d, eqn1V=%2d, eqn2V=%s*%s\n",id5b,a2,b2,g2a,m2,x,y,c,eqn1V,content(eqn2aV),factor(eqn2aV));
					);
					if(length(x5bList)<10,
						listput(x5bList,a2);
						listput(y5bList,c*g2a-2*c);
					);
				);
			);
			eqn2bV=m2*(z*z-g2a*(x*x+y*y))*(z*z-g2a*(x*x+y*y))-4*g2a*g2a+16;
			if(!polisirreducible(eqn2bV),
				if(eqn1V==2,
					\\if(has_linear_factor(eqn2bV) || dbg!=0,
					if(dbg!=0,
						printf("(%s): a=%7d, b=%7d, g=%10d, m=%14d, x=%2d, y=%2d, c=%6d, eqn1V=%2d, eqn2V=%s*%s\n",id6a,a2,b2,g2a,m2,x,y,c,eqn1V,content(eqn2aV),factor(eqn2bV));
					);
					if(length(x6aList)<10,
						listput(x6aList,a2);
						listput(y6aList,c*g2a-2*c);
					);
				);
				if(eqn1V==-2,
					\\if(has_linear_factor(eqn2bV) || dbg!=0,
					if(dbg!=0,
						printf("(%s): a=%7d, b=%7d, g=%10d, m=%14d, x=%2d, y=%2d, c=%6d, eqn1V=%2d, eqn2V=%s*%s\n",id6b,a2,b2,g2a,m2,x,y,c,eqn1V,content(eqn2aV),factor(eqn2bV));
					);
					if(length(x6bList)<10,
						listput(x6bList,a2);
						listput(y6bList,c*g2a-2*c);
					);
				);
			);
		);
		if(issquare(v2-4),
			g2b=sqrtint(v2-4);
			eqn1V=b2*(x^2-y^2)+2*a2*x*y;
			if(abs(eqn1V)!=2,
				printf("BAD: (G2B), a=%7d, b=%7d, g=%10d, m=%14d, x=%2d, y=%2d, c=%6d, eqn1V=%2d\n",a2,b2,g2b,m2,x,y,c,eqn1V);
				return();
			);
			eqn2aV=m2*(z*z-g2b*(x*x+y*y))*(z*z-g2b*(x*x+y*y))-4*g2b*g2b-16;
			if(!polisirreducible(eqn2aV),
				if(eqn1V==2,
					\\if(has_linear_factor(eqn2aV) || dbg!=0,
					if(dbg!=0,
						printf("(%s): a=%7d, b=%7d, g=%10d, m=%14d, x=%2d, y=%2d, c=%6d, eqn1V=%2d, eqn2V=%s*%s\n",id7a,a2,b2,g2b,m2,x,y,c,eqn1V,content(eqn2aV),factor(eqn2aV));
					);
					if(length(x7aList)<10,
						listput(x7aList,a2);
						listput(y7aList,c*g2b-2*c);
					);
				);
				if(eqn1V==-2,
					\\if(has_linear_factor(eqn2aV) || dbg!=0,
					if(dbg!=0,
						printf("(%s): a=%7d, b=%7d, g=%10d, m=%14d, x=%2d, y=%2d, c=%6d, eqn1V=%2d, eqn2V=%s*%s\n",id7b,a2,b2,g2b,m2,x,y,c,eqn1V,content(eqn2aV),factor(eqn2aV));
					);
					if(length(x7bList)<10,
						listput(x7bList,a2);
						listput(y7bList,c*g2b-2*c);
					);
				);
			);
			eqn2bV=m2*(z*z-g2b*(x*x+y*y))*(z*z-g2b*(x*x+y*y))-4*g2b*g2b+16;
			if(!polisirreducible(eqn2bV),
				if(eqn1V==2,
					\\if(has_linear_factor(eqn2bV) || dbg!=0,
					if(dbg!=0,
						printf("(%s): a=%7d, b=%7d, g=%10d, m=%14d, x=%2d, y=%2d, c=%6d, eqn1V=%2d, eqn2V=%s*%s\n",id8a,a2,b2,g2b,m2,x,y,c,eqn1V,content(eqn2bV),factor(eqn2bV));
					);
					if(length(x8aList)<10,
						listput(x8aList,a2);
						listput(y8aList,c*g2b-2*c);
					);
				);
				if(eqn1V==-2,
					\\if(has_linear_factor(eqn2bV) || dbg!=0,
					if(dbg!=0,
						printf("(%s): a=%7d, b=%7d, g=%10d, m=%14d, x=%2d, y=%2d, c=%6d, eqn1V=%2d, eqn2V=%s*%s\n",id8b,a2,b2,g2b,m2,x,y,c,eqn1V,content(eqn2bV),factor(eqn2bV));
					);
					if(length(x8bList)<10,
						listput(x8bList,a2);
						listput(y8bList,c*g2b-2*c);
					);
				);
			);
		);
	);
	
	check_list(x1aList,y1aList,a11,a10,b10Coeff,x,y,id1a,dbg);
	check_list(x1bList,y1bList,a11,a10,b10Coeff,x,y,id1b,dbg);
	check_list(x2aList,y2aList,a11,a10,b10Coeff,x,y,id2a,dbg);
	check_list(x2bList,y2bList,a11,a10,b10Coeff,x,y,id2b,dbg);
	check_list(x3aList,y3aList,a11,a10,b10Coeff,x,y,id3a,dbg);
	check_list(x3bList,y3bList,a11,a10,b10Coeff,x,y,id3b,dbg);
	check_list(x4aList,y4aList,a11,a10,b10Coeff,x,y,id4a,dbg);
	check_list(x4bList,y4bList,a11,a10,b10Coeff,x,y,id4b,dbg);
	
	check_list(x5aList,y5aList,a21,a20,b20Coeff,x,y,id5a,dbg);
	check_list(x5bList,y5bList,a21,a20,b20Coeff,x,y,id5b,dbg);
	check_list(x6aList,y6aList,a21,a20,b20Coeff,x,y,id6a,dbg);
	check_list(x6bList,y6bList,a21,a20,b20Coeff,x,y,id6b,dbg);
	check_list(x7aList,y7aList,a21,a20,b20Coeff,x,y,id7a,dbg);
	check_list(x7bList,y7bList,a21,a20,b20Coeff,x,y,id7b,dbg);
	check_list(x8aList,y8aList,a21,a20,b20Coeff,x,y,id8a,dbg);
	check_list(x8bList,y8bList,a21,a20,b20Coeff,x,y,id8b,dbg);
}

\\ looks like these polynomials, z2Poly, always have c^3 as the leading coefficient
\\ 27 Dec 2024 (pulled out of function above)
check_list(xList,yList,a1Coeff,a0Coeff,b0Coeff,x,y,msg,dbg=0)={
	my(aPoly,bPoly,c,eqn1V,eqn2V,g,g2Poly,gPoly,isFound,isSqr,iUB,m,mPoly,sqrCount,z,z2,z2Poly);
	
	iUB=100*1000*1000;
	if(length(xList)>9,
		z2Poly=polinterpolate(Vec(xList),Vec(yList),X);
		z2Poly=subst(z2Poly,X,a1Coeff*X1+a0Coeff);
		print("\n",msg,"(a): z2Poly=",z2Poly,", content(z2Poly)=",content(z2Poly));
		isSqr=is_poly_square(z2Poly,dbg);
		if(isSqr==1,
			aPoly=a1Coeff*X1+a0Coeff;
			bPoly=2*aPoly*x*y/(y^2-x^2)+b0Coeff;
			mPoly=aPoly*aPoly+bPoly*bPoly;
			c=x*x+y*y;
			g2Poly=mPoly*c*c-4;
			gPoly=get_sqrt(g2Poly);
			eqn1V=bPoly*(x^2-y^2)+2*aPoly*x*y;
			eqn2V=mPoly*(z2Poly-gPoly*(x^2+y^2))^2-4*gPoly^2;
			if(abs(eqn1V)==2 && abs(eqn2V)==16,
				print(msg,"(a): a=",aPoly,", b=",bPoly,", m=",mPoly,", g=",gPoly,", z^2=",z2Poly,", b(x^2-y^2)+2axy=",eqn1V,", m*(z^2-g(x^2+y^2))^2-4g^2=",eqn2V);
				sqrCount=0;
				for(i=1,iUB,
					z2=subst(z2Poly,X1,i);
					if(issquare(z2),
						sqrCount++;
						if(sqrCount<10,
							z=sqrtint(z2);
							g=subst(gPoly,X1,i);
							m=subst(mPoly,X1,i);
							print(msg,"(a): X1=",i,", g=",g,", m=",m,", z=",z);
							fullAreGAndMOK(g,m,1);
							if(fullAreGAndMOK(g,m)==1,
								isFound=simple_check_poly(subst(aPoly,X1,i),subst(bPoly,X1,i),g,x,y,z,dbg);
								if(isFound==0,
									print(msg,"(a): no value of t found for X1=",i,", g=",g,", m=",m,", z=",z);
								);
							);
						);
					);
				);
			);
			if(abs(eqn1V)!=2 || abs(eqn2V)!=16,
				print("BAD!! ",msg,"(a): a=",aPoly,", b=",bPoly,", m=",mPoly,", g=",gPoly,", z^2=",z2Poly,", b(x^2-y^2)+2axy=",eqn1V,", m*(z^2-g(x^2+y^2))^2-4g^2=",eqn2V);
			);
		);
		
		z2Poly=z2Poly+4*c;
		print("\n",msg,"(b): z2Poly=",z2Poly,", content(z2Poly)=",content(z2Poly));
		isSqr=is_poly_square(z2Poly,dbg);
		if(isSqr==1,
			eqn2V=mPoly*(z2Poly-gPoly*(x^2+y^2))^2-4*gPoly^2;
			if(abs(eqn1V)==2 && abs(eqn2V)==16,
				print(msg,"(b): a=",aPoly,", b=",bPoly,", m=",mPoly,", g=",gPoly,", z^2=",z2Poly,", b(x^2-y^2)+2axy=",eqn1V,", m*(z^2-g(x^2+y^2))^2-4g^2=",eqn2V);
				sqrCount=0;
				for(i=1,iUB,
					z2=subst(z2Poly,X1,i);
					if(issquare(z2),
						sqrCount++;
						if(sqrCount<10,
							z=sqrtint(z2);
							g=subst(gPoly,X1,i);
							m=subst(mPoly,X1,i);
							print(msg,"(b): X1=",i,", g=",g,", m=",m,", z=",z);
							fullAreGAndMOK(g,m,1);
							if(fullAreGAndMOK(g,m)==1,
								isFound=simple_check_poly(subst(aPoly,X1,i),subst(bPoly,X1,i),g,x,y,z,dbg);
								if(isFound==0,
									print(msg,"(b): no value of t found for X1=",i,", g=",g,", m=",m,", z=",z);
								);
							);
						);
					);
				);
			);
			if(abs(eqn1V)!=2 || abs(eqn2V)!=16,
				print("BAD!! ",msg,"(b): a=",aPoly,", b=",bPoly,", m=",mPoly,", g=",gPoly,", z^2=",z2Poly,", b(x^2-y^2)+2axy=",eqn1V,", m*(z^2-g(x^2+y^2))^2-4g^2=",eqn2V);
			);
		);
	);
}

\\ 25 Dec 2024
has_linear_factor(f)={
	my(fSet,fSize,p);
	
	fSet=factor(f);
	fSize=matsize(fSet)[1];
	for(i=1,fSize,
		p=fSet[i,1];
		if(poldegree(p)==1,return(1));
	);
	return(0);
}

\\ 26 Dec 2024
get_sqrt(g2Poly)={
	my(c,e,f,fSet,fSize,p);
	
	c=content(g2Poly);
	if(!issquare(c),
		print("bad content=",c," for g2Poly=",g2Poly);
		return();
	);
	fSet=factor(g2Poly);
	fSize=matsize(fSet)[1];
	f=sqrtint(c);
	for(i=1,fSize,
		p=fSet[i,1];
		e=fSet[i,2];
		if(e%2==1,
			print("non-square factor, ",p,", e=",e,", for g2Poly=",g2Poly);
			return();
		);
		f=f*p^(e/2);
	);
	return(f);
}

\\ 28 Dec 2024
is_poly_square(f,dbg=0)={
	my(c,cSqr,f1,fSet,fSize,isOK,k,m,p,r,r1);
	
	if(poldegree(f)!=1,
		print("ERROR: f=",f," must be of degree 1");
		return();
	);
	c=content(f);
	cSqr=c/core(c);
	if(dbg!=0 && cSqr!=c,
		print("is_poly_square(",f,") may not be conclusive as c=",c," is not a square");
	);
	f1=f/cSqr;
	fSet=factor(pollead(f1));
	fSize=matsize(fSet)[1];
	m=1;
	for(i=1,fSize,
		p=fSet[i,1];
		if(p>2,
			m=m*p;
		);
	);
	r=pollead(Pol(Vec(f1)%m));
	if(dbg!=0,
		print("f1=",f1,", r=",r,", type(r)=",type(r),", m=",m);
	);
	fSet=factor(abs(m));
	fSize=matsize(fSet)[1];
	for(i=1,fSize,
		p=fSet[i,1];
		if(p>2,
			k=kronecker(r,p);
			if(k==-1,
				print("f=",f,"=",cSqr,"*(",r," mod ",m,") cannot be a square mod ",p);
				return(0);
			);
			if(r%p==0 && pollead(f1)%(p*p)==0,
				r1=pollead(Pol(Vec(f1)%(p*p)));
				isOK=0;
				for(j=1,p*p,
					if((j*j)%(p*p)==r1,
						isOK=1;
					);
				);
				if(isOK==0,
					print("f=",f,"=",cSqr,"*(",r1," mod ",(p*p),") cannot be a square mod ",p*p,"=",p,"^2");
					return(0);
				);
			);
		);
	);
	return(1);
}

\\
\\ OLD functions. Initial functions before generalising them
\\

\\ 25 Dec 2024
OLD_x1_y2_search(dbg=0)={
	my(bnd,x,y);
	
	x=1;
	y=2;
	for(v=1,100000,
		if(v%100==0,print("v=",v));
		\\ from -2 on RHS in first equation
		a1=3*v+2;
		m1=25*v*v+28*v+8;
		\\ from 2 on RHS in first equation
		a2=3*v+1;
		m2=25*v*v+22*v+5;
		
		b=4*v+2; \\ regardless if from a1 or a2
		for(c=1,1000000,
			v1=m1*c*c;
			if(issquare(v1+4),
				g1a=sqrtint(v1+4);
				t1=b*(x^2-y^2)+2*a1*x*y;
				p1=m1*(z*z-g1a*(x*x+y*y))*(z*z-g1a*(x*x+y*y))-4*g1a*g1a-16;
				p2=m1*(z*z-g1a*(x*x+y*y))*(z*z-g1a*(x*x+y*y))-4*g1a*g1a+16;
				if(!polisirreducible(p1) && has_linear_factor(p1),
					printf("(G1A-P1): a=%7d, b=%7d, g=%10d, m=%14d, x=%2d, y=%2d, c=%6d, t1=%2d, p1=%s\n",a1,b,g1a,m1,x,y,c,t1,factor(p1));
				);
				if(!polisirreducible(p2) && has_linear_factor(p2),
					printf("(G1A-P2): a=%7d, b=%7d, g=%10d, m=%14d, x=%2d, y=%2d, c=%6d, t1=%2d, p2=%s\n",a1,b,g1a,m1,x,y,c,t1,factor(p2));
				);
			);
			if(issquare(v1-4),
				g1b=sqrtint(v1-4);
				t1=b*(x^2-y^2)+2*a1*x*y;
				p1=m1*(z*z-g1b*(x*x+y*y))*(z*z-g1b*(x*x+y*y))-4*g1b*g1b-16;
				p2=m1*(z*z-g1b*(x*x+y*y))*(z*z-g1b*(x*x+y*y))-4*g1b*g1b+16;
				if(!polisirreducible(p1) && has_linear_factor(p1),
					printf("(G1B-P1): a=%7d, b=%7d, g=%10d, m=%14d, x=%2d, y=%2d, c=%6d, t1=%2d, p1=%s\n",a1,b,g1b,m1,x,y,c,t1,factor(p1));
				);
				if(!polisirreducible(p2) && has_linear_factor(p2),
					printf("(G1B-P2): a=%7d, b=%7d, g=%10d, m=%14d, x=%2d, y=%2d, c=%6d, t1=%2d, p2=%s\n",a1,b,g1b,m1,x,y,c,t1,factor(p2));
				);
			);
			v2=m2*c*c;
			if(issquare(v2+4),
				g2a=sqrtint(v2+4);
				t1=b*(x^2-y^2)+2*a2*x*y;
				p1=m2*(z*z-g2a*(x*x+y*y))*(z*z-g2a*(x*x+y*y))-4*g2a*g2a-16;
				p2=m2*(z*z-g2a*(x*x+y*y))*(z*z-g2a*(x*x+y*y))-4*g2a*g2a+16;
				if(!polisirreducible(p1) && has_linear_factor(p1),
					printf("(G2A-P1): a=%7d, b=%7d, g=%10d, m=%14d, x=%2d, y=%2d, c=%6d, t1=%2d, p1=%s\n",a2,b,g2a,m2,x,y,c,t1,factor(p1));
				);
				if(!polisirreducible(p2) && has_linear_factor(p2),
					printf("(G2A-P2): a=%7d, b=%7d, g=%10d, m=%14d, x=%2d, y=%2d, c=%6d, t1=%2d, p2=%s\n",a2,b,g2a,m2,x,y,c,t1,factor(p2));
				);
			);
			if(issquare(v2-4),
				g2b=sqrtint(v2-4);
				t1=b*(x^2-y^2)+2*a2*x*y;
				p1=m2*(z*z-g2b*(x*x+y*y))*(z*z-g2b*(x*x+y*y))-4*g2b*g2b-16;
				p2=m2*(z*z-g2b*(x*x+y*y))*(z*z-g2b*(x*x+y*y))-4*g2b*g2b+16;
				if(!polisirreducible(p1) && has_linear_factor(p1),
					printf("(G2B-P1): a=%7d, b=%7d, g=%10d, m=%14d, x=%2d, y=%2d, c=%6d, t1=%2d, p1=%s\n",a2,b,g2b,m2,x,y,c,t1,factor(p1));
				);
				if(!polisirreducible(p2) && has_linear_factor(p2),
					printf("(G2B-P2): a=%7d, b=%7d, g=%10d, m=%14d, x=%2d, y=%2d, c=%6d, t1=%2d, p2=%s\n",a2,b,g2b,m2,x,y,c,t1,factor(p2));
				);
			);
		);
	);
}

\\ 25 Dec 2024
OLD_x3_y4_search(dbg=0)={
	my(bnd,x,y);
	
	x=3;
	y=4;
	for(v=1,100000,
		if(v%100==0,print("v=",v));
		\\ from -2 on RHS in first equation
		a1=7*v+3;
		b1=24*v+10;
		m1=625*v*v+522*v+109;
		\\ from 2 on RHS in first equation
		a2=7*v+4;
		b2=24*v+14;
		m2=625*v*v+728*v+212;
		
		for(c=6,10000000,
			v1=m1*c*c;
			if(issquare(v1+4),
				g1a=sqrtint(v1+4);
				t1=b1*(x^2-y^2)+2*a1*x*y;
				p1=m1*(z*z-g1a*(x*x+y*y))*(z*z-g1a*(x*x+y*y))-4*g1a*g1a-16;
				p2=m1*(z*z-g1a*(x*x+y*y))*(z*z-g1a*(x*x+y*y))-4*g1a*g1a+16;
				if(!polisirreducible(p1) && has_linear_factor(p1),
					printf("(G1A-P1): a=%7d, b=%7d, g=%10d, m=%14d, x=%2d, y=%2d, c=%6d, t1=%2d, p1=%s\n",a1,b1,g1a,m1,x,y,c,t1,factor(p1));
				);
				if(!polisirreducible(p2) && has_linear_factor(p2),
					printf("(G1A-P2): a=%7d, b=%7d, g=%10d, m=%14d, x=%2d, y=%2d, c=%6d, t1=%2d, p2=%s\n",a1,b1,g1a,m1,x,y,c,t1,factor(p2));
				);
			);
			if(issquare(v1-4),
				g1b=sqrtint(v1-4);
				t1=b1*(x^2-y^2)+2*a1*x*y;
				p1=m1*(z*z-g1b*(x*x+y*y))*(z*z-g1b*(x*x+y*y))-4*g1b*g1b-16;
				p2=m1*(z*z-g1b*(x*x+y*y))*(z*z-g1b*(x*x+y*y))-4*g1b*g1b+16;
				if(!polisirreducible(p1) && has_linear_factor(p1),
					printf("(G1B-P1): a=%7d, b=%7d, g=%10d, m=%14d, x=%2d, y=%2d, c=%6d, t1=%2d, p1=%s\n",a1,b1,g1b,m1,x,y,c,t1,factor(p1));
				);
				if(!polisirreducible(p2) && has_linear_factor(p2),
					printf("(G1B-P2): a=%7d, b=%7d, g=%10d, m=%14d, x=%2d, y=%2d, c=%6d, t1=%2d, p2=%s\n",a1,b1,g1b,m1,x,y,c,t1,factor(p2));
				);
			);
			v2=m2*c*c;
			if(issquare(v2+4),
				g2a=sqrtint(v2+4);
				t1=b2*(x^2-y^2)+2*a2*x*y;
				p1=m2*(z*z-g2a*(x*x+y*y))*(z*z-g2a*(x*x+y*y))-4*g2a*g2a-16;
				p2=m2*(z*z-g2a*(x*x+y*y))*(z*z-g2a*(x*x+y*y))-4*g2a*g2a+16;
				if(!polisirreducible(p1) && has_linear_factor(p1),
					printf("(G2A-P1): a=%7d, b=%7d, g=%10d, m=%14d, x=%2d, y=%2d, c=%6d, t1=%2d, p1=%s\n",a2,b2,g2a,m2,x,y,c,t1,factor(p1));
				);
				if(!polisirreducible(p2) && has_linear_factor(p2),
					printf("(G2A-P2): a=%7d, b=%7d, g=%10d, m=%14d, x=%2d, y=%2d, c=%6d, t1=%2d, p2=%s\n",a2,b2,g2a,m2,x,y,c,t1,factor(p2));
				);
			);
			if(issquare(v2-4),
				g2b=sqrtint(v2-4);
				t1=b2*(x^2-y^2)+2*a2*x*y;
				p1=m2*(z*z-g2b*(x*x+y*y))*(z*z-g2b*(x*x+y*y))-4*g2b*g2b-16;
				p2=m2*(z*z-g2b*(x*x+y*y))*(z*z-g2b*(x*x+y*y))-4*g2b*g2b+16;
				if(!polisirreducible(p1) && has_linear_factor(p1),
					printf("(G2B-P1): a=%7d, b=%7d, g=%10d, m=%14d, x=%2d, y=%2d, c=%6d, t1=%2d, p1=%s\n",a2,b2,g2b,m2,x,y,c,t1,factor(p1));
				);
				if(!polisirreducible(p2) && has_linear_factor(p2),
					printf("(G2B-P2): a=%7d, b=%7d, g=%10d, m=%14d, x=%2d, y=%2d, c=%6d, t1=%2d, p2=%s\n",a2,b2,g2b,m2,x,y,c,t1,factor(p2));
				);
			);
		);
	);
}

\\ (15:42) gp > \r C:\MyStuff\TexFiles\mypapers\monogenic\quartic-cyclic-monogenic\pari\xy-search.gp
\\ (15:44) gp > is_poly_square(603351125*X1 + 72804355,1)
\\ is_poly_square(603351125*X1 + 72804355) may not be conclusive as c=845 is not a square
\\ f1=3570125*X1 + 430795, r=40, type(r)=t_INT, m=65
\\ %248 = 1
\\
\\ (15:44) gp > \r C:\MyStuff\TexFiles\mypapers\monogenic\quartic-cyclic-monogenic\pari\xy-search.gp
\\ (15:46) gp > is_poly_square(603351125*X1 + 72804355,1)
\\ is_poly_square(603351125*X1 + 72804355) may not be conclusive as c=845 is not a square
\\ f1=3570125*X1 + 430795, r=40, type(r)=t_INT, m=65
\\ f=603351125*X1 + 72804355=169*(20 mod 25) cannot be a square mod 25=5^2
\\ %258 = 0
\\ the fix was to replace
\\ if(r==0 && pollead(f1)%(p*p)==0,
\\ with
\\ if(r%p==0 && pollead(f1)%(p*p)==0,
\\ 28 Dec 2024
test_sqr1(dbg=0)={
	my(actV,expV,f);
	
	f=603351125*X1 + 72804355;
	actV=is_poly_square(f,dbg);
	expV=0;
	if(actV!=expV,
		print("ERROR: actV=",actV,", expV=",expV);
	);
}

\\ 31 Dec 2024
cong_search()={
	forstep(x=213,1339,2,
		m1=x*x+1;
		v1=10*x^7+42*x^5+70*x^3+70*x-32;
		m2=m1*m1;
		v2=(10*x^7+42*x^5+70*x^3+70*x-32)%m2;
		m3=m1*m2;
		v3=(10*x^7+42*x^5+70*x^3+70*x-32)%m3;
		m4=m1*m3;
		v4=(10*x^7+42*x^5+70*x^3+70*x-32)%m4;
		for(i1=0,m1-1,
			r1=(i1*i1-v1)%m1;
			if(r1==0,
				for(i2=0,m2/m1-1,
					r2=(m2/m1)*i2+i1;
					r2Sqr=(r2*r2-v2)%m2;
					if(r2Sqr==0,
						for(i3=0,m3/m2-1,
							r3=(m3/m1)*i3+(m2/m1)*i2+i1;
							r3Sqr=(r3*r3-v3)%m3;
							if(r3Sqr==0,
								for(i4=0,m4/m3-1,
									r4=(m4/m1)*i4+(m3/m1)*i3+(m2/m1)*i2+i1;
									r4Sqr=(r4*r4-v4)%m4;
									if(r4Sqr==0,
										printf("found! x=%4d, r4=%25d, r4^2 mod m4=%25d, v4 mod m4=%25d, r4/x^8=%12.8f\n",x,r4,(r4*r4)%m4,v4,1.0*r4/x^8);
									);
								);
							);
						);
					);
				);
			);
		);
		if(x<12,
			for(i1=1,m4,
				i1Sqr=(i1*i1-v4)%m4;
				if(i1Sqr==0,
					printf("brute force: x=%4d, r4=%12d, r4^2 mod m4=%12d, v4 mod m4=%12d\n",x,i1,(i1*i1)%m4,v4);
				);
			);
		);
	);
}