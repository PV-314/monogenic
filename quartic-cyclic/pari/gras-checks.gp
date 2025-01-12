read("gras-utils.gp");
\\ Gras Besancon 1981 paper Theoreme 1
\\ main entry point to other functions
\\ check_family1() at the end may also be of interest
\\ 17 May 2024
thm1_check(dbg=0)={
	my(aUB,bUB,isOK,m,v1,v1s,v2,xLB,xUB,y1,y2);

	aUB=1000; \\ 1000*1000;
	bUB=100*1000;
	xLB=0;
	xUB=200;
	for(a=2,aUB,
		if(a%10==0,print("a=",a));
		forstep(b=2,bUB,2,
			m=a*a+b*b;
			isOK=isMOK(m);
			if(isOK,
				if(dbg!=0 && m==200,
					print("m=",m,", isOK=",isOK);
				);
				for(x=xLB,xUB,
					\\ from f1=b*(x*x-y*y)+2*a*x*y-(+2);
					v1=a*a*x*x+b*b*x*x-2*b;
					if(issquare(v1),
						v1s=sqrtint(v1);
						y1=(2*a*x+2*v1s)/2/b;
						if(denominator(y1)==1,
							check_eqn2_by_g(a,b,m,x,y1,dbg);
							check_eqn2_by_z(a,b,m,x,y1,dbg);
						);
						y2=(2*a*x-2*v1s)/2/b;
						if(denominator(y2)==1,
							\\check_eqn2_by_g(a,b,m,x,y2,dbg);
							check_eqn2_by_z(a,b,m,x,y2,dbg);
						);
					);

					\\ from f2=b*(x*x-y*y)+2*a*x*y-(-2);
					v2=a*a*x*x+b*b*x*x+2*b;
					if(issquare(v2),
						v2s=sqrtint(v2);
						y1=(2*a*x+2*v2s)/2/b;
						if(denominator(y1)==1,
							check_eqn2_by_g(a,b,m,x,y1,dbg);
							check_eqn2_by_z(a,b,m,x,y1,dbg);
						);
						y2=(2*a*x-2*v2s)/2/b;
						if(denominator(y2)==1,
							check_eqn2_by_g(a,b,m,x,y2,dbg);
							check_eqn2_by_z(a,b,m,x,y2,dbg);
						);
					);
				);
			);
		);
	);
}

\\ used to find families for other pairs, (x,y)
thm1_check_x_y_all(aUB=10000,dbg=0)={
	for(x=1,30,
	for(y=x+1,x+20,
		if(x!=y && gcd(x,y)==1,
			thm1_check_for_x_y(x,y,aUB,dbg);
		);
	);
	);
}

\\ search for solutions when x=1 and y=2
\\ no further examples besides a=1, b=2,... and
\\ FOUND        : for a=   5, b=   6, g= 1523, m=      61, x=   1, y=   2, z=  85
\\ found
\\ 20 May 2024
thm1_check_for_x_y(x,y,aUB=1000,dbg=0)={
	my(b,b1,b2,isOK,m,x2MinusY2);
	
	print("using aUB=",aUB);
	x2MinusY2=x*x-y*y;
	for(a=0,aUB,
		if(a%1000==0,
			print("for x=",x,", y=",y,": starting a=",a);
		);
		\\ from first equation in equation (1) in Gras, b=(\pm 2 - 2*a*x*y)/(x^2-y^2)
		\\ apply that here. Recall that b must be even too (hence mod 6=2*(2^2-1^2))
		b1=2-2*a*x*y;
		if(b1%x2MinusY2==0,
			b=b1/x2MinusY2;
			if(b%2==0,
				m=a*a+b*b;
				isOK=isMOK(m);
				if(isOK,
					check_eqn2_by_g(a,b,m,x,y);
					check_eqn2_by_z(a,b,m,x,y);
				);
			);
		);
		b2=-2-2*a*x*y;
		if(b2%x2MinusY2==0,
			b=b2/x2MinusY2;
			if(b%2==0,
				m=a*a+b*b;
				isOK=isMOK(m);
				if(isOK,
					check_eqn2_by_g(a,b,m,x,y);
					check_eqn2_by_z(a,b,m,x,y);
				);
			);
		);
	);
}

\\ looking for real fields, so chi(-1)=+1
\\ eqn2 is the second equation in equation (1) in Gras' paper (part of her Theoreme 1)
check_eqn2_by_g(a,b,m,x,y,dbg=0)={
	my(eqn1Value,eqn2Value,gLB,gUB,isOK,v1,v1Sqr,x2,y2,z,zSqr);
	
	if(dbg!=0,
		print("check_eqn2_by_g(): checking a=",a,", b=",b,", m=",m,", x=",x,", y=",y);
	);
	gLB=1;
	gUB=1000*1000;
	x2=x*x;
	y2=y*y;
	for(g=gLB,gUB,
		if(dbg!=0,
			if(a==41,
				print("check_eqn2_by_g(): a=",a,", b=",b,", g=",g,", x=",x,", y=",y);
			);
		);

		\\ these conditions come from the end of the first paragraph of Section 1 in Gras' paper
		isOK=areGAndMOK(g,m);
		if(isOK,
			if(dbg!=0,
				if(a==41,
					print("a=",a,", b=",b,", g=",g,", x=",x,", y=",y,", m*g*g+4*m=",m*g*g+4*m);
				);
			);
			\\ from m*(z^2-g*(x^2+y^2))^2-4*g*g-(+16)
			v1Sqr=4*g*g+16;
			if(v1Sqr%m==0,
				v1Sqr=v1Sqr/m;
				if (issquare(v1Sqr),
					v1=sqrtint(v1Sqr);
					zSqr=v1+g*(x2+y2);
					if (issquare(zSqr),
						z=sqrtint(zSqr);
						eqn1Value=b*(x2-y2)+2*a*x*y;
						eqn2Value=m*(z*z-g*(x2+y2))*(z*z-g*(x2+y2))-4*g*g;
						printf("(G:%2d:%3d:a1): a=%5d, b=%5d, m=%10d, g=%7d, x=%5d, y=%4d, z=%5d\n",eqn1Value,eqn2Value,a,b,m,g,x,y,z);
						check_poly(a,b,g,x,y,z);
					);
					\\ the other square root of v1Sqr
					v1=-v1;
					zSqr=v1+g*(x2+y2);
					if (issquare(zSqr),
						z=sqrtint(zSqr);
						eqn1Value=b*(x2-y2)+2*a*x*y;
						eqn2Value=m*(z*z-g*(x2+y2))*(z*z-g*(x2+y2))-4*g*g;
						printf("(G:%2d:%3d:a2): a=%5d, b=%5d, m=%10d, g=%7d, x=%5d, y=%4d, z=%5d\n",eqn1Value,eqn2Value,a,b,m,g,x,y,z);
						check_poly(a,b,g,x,y,z);
					);
				);
			);
			v2Sqr=4*g*g-16;
			if(v2Sqr%m==0,
				v2Sqr=v2Sqr/m;
				if (issquare(v2Sqr),
					v2=sqrtint(v2Sqr);
					z2=v2+g*(x2+y2);
					if (issquare(z2),
						z=sqrtint(z2);
						eqn1Value=b*(x2-y2)+2*a*x*y;
						eqn2Value=m*(z*z-g*(x2+y2))*(z*z-g*(x2+y2))-4*g*g;
						printf("(G:%2d:%3d:b1): a=%5d, b=%5d, m=%10d, g=%7d, x=%5d, y=%4d, z=%5d\n",eqn1Value,eqn2Value,a,b,m,g,x,y,z);
						check_poly(a,b,g,x,y,z);
					);

					v2=-v2;
					z2=v2+g*(x2+y2);
					if (issquare(z2),
						z=sqrtint(z2);
						eqn1Value=b*(x2-y2)+2*a*x*y;
						eqn2Value=m*(z*z-g*(x2+y2))*(z*z-g*(x2+y2))-4*g*g;
						printf("(G:%2d:%3d:b2): a=%5d, b=%5d, m=%10d, g=%7d, x=%5d, y=%4d, z=%5d\n",eqn1Value,eqn2Value,a,b,m,g,x,y,z);
						check_poly(a,b,g,x,y,z);
					);
				);
			);
		);
	);
}

\\ looking for real fields, so chi(-1)=+1
\\ eqn2 is the second equation in equation (1) in Gras' paper (part of her Theoreme 1)
check_eqn2_by_z(a,b,m,x,y,dbg=0)={
	my(gLB,gUB,isOK,v1,v2,v2Sqr,v3,v4,v4Sqr,v5);
	
	x2=x*x;
	y2=y*y;
	zLB=1;
	zUB=5*1000*1000;
	tmp1a=4*m*x2*x2+8*m*x2*y2+4*m*y2*y2-16;
	for(z=zLB,zUB,
		z2=z*z;
		tmp2=tmp1a+m*z2*z2;
		if (issquare(tmp2),
			disc=sqrtint(tmp2);
			gNum1=2*m*x2*z2+2*m*y2*z2+4*disc;
			gDen=2*(m*x2*x2+2*m*x2*y2+m*y2*y2-4);
			if(gDen!=0 && gNum1%gDen==0,
				g=gNum1/gDen;
				isOK=areGAndMOK(g,m);
				if(isOK,
					eqn1Value=b*(x2-y2)+2*a*x*y;
					eqn2Value=m*(z2-g*(x2+y2))*(z2-g*(x2+y2))-4*g*g;
					printf("(Z:%2d:%3d:a1): a=%5d, b=%5d, m=%9d, g=%7d, x=%5d, y=%4d, z=%5d\n",eqn1Value,eqn2Value,a,b,m,g,x,y,z);
					check_poly(a,b,g,x,y,z);
				);
			);
			gNum2=2*m*x2*z2+2*m*y2*z2-4*disc;
			if (gDen!=0 && gNum2%gDen==0,
				g=gNum2/gDen;
				isOK=areGAndMOK(g,m);
				if(isOK,
					eqn1Value=b*(x2-y2)+2*a*x*y;
					eqn2Value=m*(z2-g*(x2+y2))*(z2-g*(x2+y2))-4*g*g;
					printf("(Z:%2d:%3d:a2): a=%5d, b=%5d, m=%9d, g=%7d, x=%5d, y=%4d, z=%5d\n",eqn1Value,eqn2Value,a,b,m,g,x,y,z);
					check_poly(a,b,g,x,y,z);
				);
			);
		);
	);
	tmp1b=-4*m*x2*x2-8*m*x2*y2-4*m*y2*y2+16;
	for(z=zLB,zUB,
		z2=z*z;
		tmp2=tmp1b+m*z2*z2;
		if (issquare(tmp2),
			disc=sqrtint(tmp2);
			gNum1=2*m*x2*z2+2*m*y2*z2+4*disc;
			gDen=2*(m*x2*x2+2*m*x2*y2+m*y2*y2-4);
			if(gDen!=0 && gNum1%gDen==0,
				g=gNum1/gDen;
				isOK=areGAndMOK(g,m);
				if(isOK,
					eqn1Value=b*(x2-y2)+2*a*x*y;
					eqn2Value=m*(z*z-g*(x2+y2))*(z*z-g*(x2+y2))-4*g*g;
					printf("(Z:%2d:%3d:b1): a=%5d, b=%5d, m=%9d, g=%7d, x=%5d, y=%4d, z=%5d\n",eqn1Value,eqn2Value,a,b,m,g,x,y,z);
					check_poly(a,b,g,x,y,z);
				);
			);
			gNum2=2*m*x2*z2+2*m*y2*z2-4*disc;
			if (gDen!=0 && gNum2%gDen==0,
				g=gNum2/gDen;
				isOK=areGAndMOK(g,m);
				if(isOK,
					eqn1Value=b*(x2-y2)+2*a*x*y;
					eqn2Value=m*(z*z-g*(x2+y2))*(z*z-g*(x2+y2))-4*g*g;
					printf("(Z:%2d:%3d:b2): a=%5d, b=%5d, m=%9d, g=%7d, x=%5d, y=%4d, z=%5d\n",eqn1Value,eqn2Value,a,b,m,g,x,y,z);
					check_poly(a,b,g,x,y,z);
				);
			);
		);
	);
}
