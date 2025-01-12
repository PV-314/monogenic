\\ for showing results for z=-1,0,1, where the fields are not totally real
\\ but other values of z can be tested too

\\ 4 Jan 2025
z_check(z,dbg=0)={
	my(df,f,nf,rRoots);


	f=x^4-z*x^3+(-3/8*z^6+2*z^4-37/8*z^2+4)*x^2+(19/16*z^7-9/2*z^5+135/16*z^3-6*z-1/8*z^9)*x+2+9/4*z^6-107/128*z^8+1/4*z^2+5/32*z^10-3/256*z^12-719/256*z^4;
	printf("for z=%1d, f=%s\n",z,f);
	rRoots=polrootsreal(f);
	if(length(rRoots)==poldegree(f),
		printf("for z=%1d, polroots(f)=%10.8f\n",z,rRoots);
	);
	if(length(rRoots)!=poldegree(f),
		printf("for z=%1d, polroots(f)=%10.8f\n",z,polroots(f));
	);
	print("for z=",z,", polisirreducible(f)=",polisirreducible(f));
	if(!polisirreducible(f),
		print("ERROR: f is not irreducible for z=",z,": f=",factor(f));
		return();
	);
	print("for z=",z,", polgalois(f)=",polgalois(f));
	df=poldisc(f);
	nf=bnfinit(f,1);
	if(nf.disc!=df,
		print("for z=",z,", f is not monogenic. disc(f)/nf.disc=",df/nf.disc);
	);
	if(nf.disc==df,
		print("for z=",z,", f is monogenic. disc(f)=",df);
	);
}