\\ 12 Jan 2025
fz_integer_coeff_check(dbg=0)={
	my(f);

	f=x^4-z*x^3+(-3/8*z^6+2*z^4-37/8*z^2+4)*x^2+(19/16*z^7-9/2*z^5+135/16*z^3-6*z-1/8*z^9)*x+2+9/4*z^6-107/128*z^8+1/4*z^2+5/32*z^10-3/256*z^12-719/256*z^4;
	f_integer_coeff_check(f,0);
}

\\ for checking f_z(X-1/2)
\\ 12 Jan 2025
fz12_integer_coeff_check(dbg=0)={
	my(f);

	f=x^4-z*x^3+(-3/8*z^6+2*z^4-37/8*z^2+4)*x^2+(19/16*z^7-9/2*z^5+135/16*z^3-6*z-1/8*z^9)*x+2+9/4*z^6-107/128*z^8+1/4*z^2+5/32*z^10-3/256*z^12-719/256*z^4;
	f=subst(f,x,x-1/2);
	print(f);
	f_integer_coeff_check(f,1);
}

\\ aOrB is 0 if for part (a) and 1 for part (b)
\\ 12 Jan 2025
f_integer_coeff_check(f,aOrB,dbg=0)={
	my(cf0,cf0a,cf0Denom,cf0Numer,cf1,cf1a,cf1Denom,cf1Numer,cf2,cf2a,cf2Denom,cf2Numer);

	if(aOrB!=0 && aOrB!=1,
		print("ERROR: aOrB must be 0 or 1, but is",aOrB);
		return();
	);
	cf2=polcoef(f,2,x);
	cf2Denom=denominator(content(cf2));
	cf2Numer=cf2Denom*cf2;
	if(aOrB==0,
		cf2a=Pol(Vec(subst(cf2Numer,z,2*z1))%cf2Denom,z1);
		print("for z even=  2*z1, numer of cf2 mod   ",cf2Denom,"(the denominator) is ",cf2a);
		cf2a=Pol(Vec(subst(cf2Numer,z,2*z1+1))%cf2Denom,z1);
		print("for z odd =2*z1+1, numer of cf2 mod   ",cf2Denom,"(the denominator) is ",cf2a);
	);

	if(aOrB==1,
		cf2a=Pol(Vec(subst(cf2Numer,z,4*z1))%cf2Denom,z1);
		print("for z even=  4*z1, numer of cf2 mod   ",cf2Denom,"(the denominator) is ",cf2a);
		cf2a=Pol(Vec(subst(cf2Numer,z,4*z1+2))%cf2Denom,z1);
		print("for z even=4*z1+2, numer of cf2 mod   ",cf2Denom,"(the denominator) is ",cf2a);
		cf2a=Pol(Vec(subst(cf2Numer,z,2*z1+1))%cf2Denom,z1);
		print("for z odd =2*z1+1, numer of cf2 mod   ",cf2Denom,"(the denominator) is ",cf2a);
	);

	cf1=polcoef(f,1,x);
	cf1Denom=denominator(content(cf1));
	cf1Numer=cf1Denom*cf1;
	if(aOrB==0,
		cf1a=Pol(Vec(subst(cf1Numer,z,2*z1))%cf1Denom,z1);
		print("for z even=  2*z1, numer of cf1 mod  ",cf1Denom,"(the denominator) is ",cf1a);
		cf1a=Pol(Vec(subst(cf1Numer,z,2*z1+1))%cf1Denom,z1);
		print("for z odd =2*z1+1, numer of cf1 mod  ",cf1Denom,"(the denominator) is ",cf1a);
	);
	if(aOrB==1,
		cf1a=Pol(Vec(subst(cf1Numer,z,4*z1))%cf1Denom,z1);
		print("for z even=  4*z1, numer of cf1 mod  ",cf1Denom,"(the denominator) is ",cf1a);
		cf1a=Pol(Vec(subst(cf1Numer,z,4*z1+2))%cf1Denom,z1);
		print("for z even=4*z1+2, numer of cf1 mod  ",cf1Denom,"(the denominator) is ",cf1a);
		cf1a=Pol(Vec(subst(cf1Numer,z,2*z1+1))%cf1Denom,z1);
		print("for z odd =2*z1+1, numer of cf1 mod  ",cf1Denom,"(the denominator) is ",cf1a);
	);

	cf0=polcoef(f,0,x);
	cf0Denom=denominator(content(cf0));
	cf0Numer=cf0Denom*cf0;
	if(aOrB==0,
		cf0a=Pol(Vec(subst(cf0Numer,z,2*z1))%cf0Denom,z1);
		print("for z even=  2*z1, numer of cf0 mod ",cf0Denom,"(the denominator) is ",cf0a);
		cf0a=Pol(Vec(subst(cf0Numer,z,2*z1+1))%cf0Denom,z1);
		print("for z odd =2*z1+1, numer of cf0 mod ",cf0Denom,"(the denominator) is ",cf0a);
	);
	if(aOrB==1,
		cf0a=Pol(Vec(subst(cf0Numer,z,4*z1))%cf0Denom,z1);
		print("for z even=  4*z1, numer of cf0 mod ",cf0Denom,"(the denominator) is ",cf0a);
		cf0a=Pol(Vec(subst(cf0Numer,z,4*z1+2))%cf0Denom,z1);
		print("for z even=4*z1+2, numer of cf0 mod ",cf0Denom,"(the denominator) is ",cf0a);
		cf0a=Pol(Vec(subst(cf0Numer,z,2*z1+1))%cf0Denom,z1);
		print("for z odd =2*z1+1, numer of cf0 mod ",cf0Denom,"(the denominator) is ",cf0a);
	);
}

\\ 12 Jan 2025
galois_psi_check(dbg=0)={
	my(a,b,c,d,discR2,f,g1,g2,r,r2,rFact,s);
	
	f=x^4-(z^2-2)*(z^4-4*z^2+8)*x^2+(z^2-2)^2*(z^4-4*z^2+8);
	a=polcoef(f,3,x);
	b=polcoef(f,2,x);
	c=polcoef(f,1,x);
	d=polcoef(f,0,x);
	r=x^3-b*x^2+(a*c-4*d)*x-(a*a*d-4*b*d+c*c);
	rFact=factor(r);
	rSize=matsize(rFact)[1];
	print("r=",rFact);
	if(rSize!=2,
		print("r should have exactly two factors, but r=",rFact);
		return();
	);
	r1=rFact[1,1];
	r2=rFact[2,1];
	if(poldegree(r1,x)!=1,
		r1=rFact[2,1];
		r2=r1;
	);
	if((poldegree(r1,x)!=1 && poldegree(r2,x)!=1) || (poldegree(r1,x)!=2 && poldegree(r2,x)!=2),
		print("r should have exactly two factors, one linear, the other quadratic, but r=",rFact);
		return();
	);
	s=-polcoef(r1,0,x)/polcoef(r1,1,x);
	discR2=poldisc(r2,x);
	print("disc(quadratic factor)=",content(discR2),"*",factor(discR2));
	
	print("\ns=",s,", d=",d);
	g1=x^2-s*x+d;
	g2=x^2+a*x+b-s;
	print("g1=x^2-s*x+d=",g1);
	print("disc(g1)=",factor(poldisc(g1,x)));
	print("g2=",g2,", disc(g2)=",factor(poldisc(g2,x)));
}
