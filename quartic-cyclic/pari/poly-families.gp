read("gras-utils.gp");

\\ the family in main theorem
\\ 18 May 2024
check_fz_polys(dbg=0)={
	my(a,av,b,bv,c,dc,f,fv,g,gv,m,mv,t,x,y);

	\\ from my family:
	kill(z);
	a=z*z-2;
	b=2;
	g=a;
	t=z;
	x=1;
	y=0;
	f=calc_poly(a,b,g,t,x,y,z);
	print(f);
	c=content(f);
	dc=denominator(c);
	print("X^3 coeff:",dc*polcoef(f,3,X)," / ",dc,"\n");
	print("X^2 coeff:",dc*polcoef(f,2,X)," / ",dc,"\n");
	print("X^1 coeff:",dc*polcoef(f,1,X)," / ",dc,"\n");
	print("X^0 coeff:",dc*polcoef(f,0,X)," / ",dc,"\n");
	c=0;
	for(zv=1,5000,
		av=zv*zv-2;
		bv=2;
		gv=av;
		mv=av^2+bv^2;
		if(issquarefree(gv) && issquarefree(mv/gcd(4,mv)),
			fv=subst(f,z,zv);
			if(zv%4==2,
				fv=subst(subst(f,z,zv),X,X-1/2);
			);
			if(denominator(content(fv))!=1,
				print("BAD: z=",zv,", f=",fv);
				return();
			);
			if(polisirreducible(fv) && denominator(content(fv))==1,
				nf=bnfinit(fv,1);
				\\print("   poldisc(f)=",poldisc(f),", nf.disc=",nf.disc);
				if(poldisc(fv)==nf.disc,
					c=c+1;
					printf("MONOGENIC, c=%3d, z=%3d, poldisc(f)=%12d, nf.disc=%12d\n",c,zv,poldisc(fv),nf.disc);
				);
			);
		);
	);
}

\\ Example (1) in the "Further Families" section of paper
\\ 24 Dec 2024
check_further_family1(dbg=0)={
	my(a,b,c,dc,f,g,m,t,x,y);

	kill(z);
	a=z*z+2;
	b=2;
	g=a;
	t=z;
	x=1;
	y=0;
	f=calc_poly(a,b,g,t,x,y,z);
	print(f);
	c=content(f);
	dc=denominator(c);
	print("X^3 coeff:",dc*polcoef(f,3,X)," / ",dc,"\n");
	print("X^2 coeff:",dc*polcoef(f,2,X)," / ",dc,"\n");
	print("X^1 coeff:",dc*polcoef(f,1,X)," / ",dc,"\n");
	print("X^0 coeff:",dc*polcoef(f,0,X)," / ",dc,"\n");
	for(z=1,5,
		a=z*z+2;
		b=2;
		g=a;
		m=a^2+b^2;
		if(issquarefree(g) && issquarefree(m/gcd(4,m)),
			print("\nz=",z);
			check_poly(a,b,g,x,y,z,dbg);
		);
	);
}

\\ Example (2) in the "Further Families" section of paper
\\ 24 Dec 2024
check_further_family2(dbg=0)={
	my(a,b,c,dc,f,g,m,t,x,y,z);
	
	kill(v);
	x=1;
	y=1;
	z=4*v+2;
	a=1;
	g=8*v^2+8*v+4;
	b=g/2;
	t=z;
	f=calc_poly(a,b,g,t,x,y,z);
	print(f);
	c=content(f);
	dc=denominator(c);
	print("X^3 coeff:",dc*polcoef(f,3,X)," / ",dc,"\n");
	print("X^2 coeff:",dc*polcoef(f,2,X)," / ",dc,"\n");
	print("X^1 coeff:",dc*polcoef(f,1,X)," / ",dc,"\n");
	print("X^0 coeff:",dc*polcoef(f,0,X)," / ",dc,"\n");
	for(v=1,10,
		z=4*v+2;
		a=1;
		g=8*v^2+8*v+4;
		b=g/2;
		t=z;
		f=calc_poly(a,b,g,t,x,y,z);
		m=a^2+b^2;
		if(issquarefree(g/gcd(2^20,g)) && issquarefree(m/gcd(2^20,m)),
			print("\nz=",z);
			check_poly(a,b,g,x,y,z,dbg);
		);
	);
}

\\ 17 Dec 2024
check_x3y4_examples(dbg=0)={
	my(a,b,g,m,x,y);
	
	a=  466;b= 1598;m=   2770760;g=  41614;x=    3;y=   4;z= 1020;
	check_poly(a,b,g,x,y,z,dbg);
	a=  717;b= 2458;m=   6555853;g=  64011;x=    3;y=   4;z= 1265;
	check_poly(a,b,g,x,y,z,dbg);

	a= 1550; b= 5314; m=  30641096; g= 138386; x=    3; y=   4; z= 1860;
	check_poly(a,b,g,x,y,z,dbg);
	a= 1985; b= 6806; m=  50261861; g= 177239; x=    3; y=   4; z= 2105;
	check_poly(a,b,g,x,y,z,dbg);

	a= 7697; b=26390; m= 755675909; g= 687239; x=    3; y=   4; z= 4145;
	check_poly(a,b,g,x,y,z,dbg);
	a= 8634; b=29602; m= 950824360; g= 770886; x=    3; y=   4; z= 4390;
	check_poly(a,b,g,x,y,z,dbg);

	a=11133; b=38170; m=1580892589; g= 994011; x=    3; y=   4; z= 4985;
	check_poly(a,b,g,x,y,z,dbg);
	a=12254; b=42014; m=1915336712; g=1094114; x=    3; y=   4; z= 5230;
	check_poly(a,b,g,x,y,z,dbg);

	a=23678; b=81182; m=7151164808; g=2114114; x=    3; y=   4; z= 7270;
	check_poly(a,b,g,x,y,z,dbg);
	a=25301; b=86746; m=8165009117; g=2259011; x=    3; y=   4; z= 7515;
	check_poly(a,b,g,x,y,z,dbg);
}

\\ 18 Dec 2024
check_x3y4_family1(dbg=0)={
	my(a,b,c,dc,f,g,m,x,y,z);

	kill(t);
	kill(z1);
	x=3;
	y=4;
	z=3125*z1+1020;
	a=7/15625*z^2-62/625;
	b=24/15625*z^2-34/625;
	m=a^2+b^2;
	g=z^2/25-2;
	f=calc_poly(a,b,g,t,x,y,z);
	\\print(f);
	c=content(f);
	dc=denominator(c);
	print("X^3 coeff:",dc*polcoef(f,3,X)," / ",dc,"\n");
	print("X^2 coeff:",dc*polcoef(f,2,X)," / ",dc,"\n");
	print("X^1 coeff:",dc*polcoef(f,1,X)," / ",dc,"\n");
	print("X^0 coeff:",dc*polcoef(f,0,X)," / ",dc,"\n");
	for(z1=1,10,
		z=3125*z1+1020;
		a=7/15625*z^2-62/625;
		b=24/15625*z^2-34/625;
		m=a^2+b^2;
		g=z^2/25-2;
		if(z1%4==1,
			t=3;
			f=calc_poly(a,b,g,t,x,y,z);
		);
		if(z1%4==3,
			t=1;
			f=calc_poly(a,b,g,t,x,y,z);
		);
		if(z1%2==0,
			t=4;
			f=calc_poly(a,b,g,t,x,y,z);
		);
			
		if(issquarefree(g/gcd(2^20,g)) && issquarefree(m/gcd(2^20,m)),
			print("z1=",z1,", m=",m,", t=",t,",z=",z,", (t+z)/2=",((t+z)/2)%2,", gx=",(g*x)%2);
			print("f=",f);
			if(denominator(content(f))!=1,
				print("BAD!");
				return();
			);
		);
	);
}

\\ 18 Dec 2024
check_x3y4_family2(dbg=0)={
	my(a,b,c,dc,f,g,m,x,y,z);

	kill(t);
	kill(z1);
	x=3;
	y=4;
	z=3125*z1+1265;
	a=7/15625*z^2+62/625;
	b=24/15625*z^2+34/625;
	m=a^2+b^2;
	g=z^2/25+2;
	f=calc_poly(a,b,g,t,x,y,z);
	\\print(f);
	c=content(f);
	dc=denominator(c);
	print("X^3 coeff:",dc*polcoef(f,3,X)," / ",dc,"\n");
	print("X^2 coeff:",dc*polcoef(f,2,X)," / ",dc,"\n");
	print("X^1 coeff:",dc*polcoef(f,1,X)," / ",dc,"\n");
	print("X^0 coeff:",dc*polcoef(f,0,X)," / ",dc,"\n");
	for(z1=1,10,
		z=3125*z1+1265;
		a=7/15625*z^2+62/625;
		b=24/15625*z^2+34/625;
		m=a^2+b^2;
		g=z^2/25+2;

		if(z1%2==1,
			t=0;
			f=calc_poly(a,b,g,t,x,y,z);
		);
		if(z1%4==0,
			t=3;
			f=calc_poly(a,b,g,t,x,y,z);
		);
		if(z1%4==2,
			t=1;
			f=calc_poly(a,b,g,t,x,y,z);
		);
			
		if(issquarefree(g/gcd(2^20,g)) && issquarefree(m/gcd(2^20,m)),
			print("z1=",z1,", m=",m,", t=",t,",z=",z,", (t+z)/2=",((t+z)/2)%2,", gx=",(g*x)%2);
			print("f=",f);
			if(denominator(content(f))!=1,
				print("BAD!");
				return();
			);
		);
	);
}
	
\\ 18 Dec 2024
check_x3y4_family3(dbg=0)={
	my(a,b,c,dc,f,g,m,x,y,z);

	kill(t);
	kill(z1);
	x=3;
	y=4;
	z=3125*z1+(3125-1265);
	a=7/15625*z^2+62/625;
	b=24/15625*z^2+34/625;
	m=a^2+b^2;
	g=z^2/25+2;
	f=calc_poly(a,b,g,t,x,y,z);
	\\print(f);
	c=content(f);
	dc=denominator(c);
	print("X^3 coeff:",dc*polcoef(f,3,X)," / ",dc,"\n");
	print("X^2 coeff:",dc*polcoef(f,2,X)," / ",dc,"\n");
	print("X^1 coeff:",dc*polcoef(f,1,X)," / ",dc,"\n");
	print("X^0 coeff:",dc*polcoef(f,0,X)," / ",dc,"\n");
	
	for(z1=1,10,
		z=3125*z1+(3125-1265);
		a=7/15625*z^2+62/625;
		b=24/15625*z^2+34/625;
		m=a^2+b^2;
		g=z^2/25+2;
		if(z1%4==1,
			t=3;
			f=calc_poly(a,b,g,t,x,y,z);
		);
		if(z1%4==3,
			t=1;
			f=calc_poly(a,b,g,t,x,y,z);
		);
		if(z1%2==0,
			t=0;
			f=calc_poly(a,b,g,t,x,y,z);
		);
		if(issquarefree(g/gcd(2^20,g)) && issquarefree(m/gcd(2^20,m)),
			print("z1=",z1,", m=",m,", t=",t,",z=",z);
			print("f=",f);
			if(denominator(content(f))!=1,
				print("BAD!");
				return();
			);
		);
	);
}

\\ 18 Dec 2024
check_x3y4_family4(dbg=0)={
	my(a,b,c,dc,f,g,m,x,y,z);

	kill(t);
	kill(z1);
	x=3;
	y=4;
	z=3125*z1+3125-1020;
	a=7/15625*z^2-62/625;
	b=24/15625*z^2-34/625;
	m=a^2+b^2;
	g=z^2/25-2;
	f=calc_poly(a,b,g,t,x,y,z);
	\\print(f);
	c=content(f);
	dc=denominator(c);
	print("X^3 coeff:",dc*polcoef(f,3,X)," / ",dc,"\n");
	print("X^2 coeff:",dc*polcoef(f,2,X)," / ",dc,"\n");
	print("X^1 coeff:",dc*polcoef(f,1,X)," / ",dc,"\n");
	print("X^0 coeff:",dc*polcoef(f,0,X)," / ",dc,"\n");
	
	for(z1=1,10,
		z=3125*z1+3125-1020;
		a=7/15625*z^2-62/625;
		b=24/15625*z^2-34/625;
		m=a^2+b^2;
		g=z^2/25-2;
		if(z1%2==1,
			t=0;
			f=calc_poly(a,b,g,t,x,y,z);
		);
		if(z1%4==2,
			t=1;
			f=calc_poly(a,b,g,t,x,y,z);
		);
		if(z1%4==0,
			t=3;
			f=calc_poly(a,b,g,t,x,y,z);
		);
		if(issquarefree(g/gcd(2^20,g)) && issquarefree(m/gcd(2^20,m)),
			print("z1=",z1,", m=",m,", t=",t,",z=",z);
			print("f=",f);
			if(denominator(content(f))!=1,
				print("BAD!");
				return();
			);
			\\simple_check_poly(a,b,g,x,y,z,dbg);
		);
	);
}
