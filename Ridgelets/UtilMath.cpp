#include "pch.h"
#include "UtilMath.h"

void UtilMath::spiralsample(MatrixType& u, unsigned flg, unsigned N)
{
	//z=(1-1/N:-2/N:1/N-1)';
	MatrixType z(N, 1);
	double val = 1.0 - (1.0 / N);
	for (unsigned int i = 0; i < N; ++i, val += (-2.0 / N))
	{
		z(i) = val;
	}

	//r=sqrt(1-z.*z);
	MatrixType r(N, 1);
	r = (1.0 - z.cwiseProduct(z).array()).sqrt();

	MatrixType Long(N, 1);
	switch (flg)
	{
	case 1:
	{
		// Long=[0; cumsum((3.6/sqrt(N))./r(1:N-1))];
		double sqrtN = sqrt(N);
		Long(0) = 0.0;
		for (unsigned i = 1; i < N; ++i)
		{
			Long(i) = Long(i - 1) + ((3.6 / sqrtN) / r(i - 1));
		}
	}
	break;
	case 2:
	{
		// Long=(pi*(3-sqrt(5)))*(0:N-1)';
		double pref_const = UtilMath::PI * (3.0 - sqrt(5.0));
		Long = pref_const * VectorXd::LinSpaced(N, 0, N - 1);
	}
	break;
	default:
		cerr << "Invalid Sampling Option \n";
		throw;
	}

	//u=[r.*cos(long) r.*sin(long) z];
	//u=u./repmat(sqrt(sum(u.^2,2)),[1 3]);

	u.col(0) = r.array() * Long.array().cos();
	u.col(1) = r.array() * Long.array().sin();
	u.col(2) = z.array();
	for (unsigned i = 0; i < N; ++i)
		u.row(i).normalize();

	//return u;
}

void UtilMath::fura(MatrixType& Lmd, unsigned n)
{
	//from fura.m

	//Lmd=ones(n+1,1);

	//Lmd(2:2:n+1)=0;
	Lmd(0) = 1.0;
	double c;

	//for k=2:2:n
	for (unsigned k = 2; k < n + 1; k += 2) {
		//Lmd(k+1)=-Lmd(k-1)*(k-1)/k;
		c = (double)k;
		Lmd.row(k) = -((c - 1.0) / c)*Lmd.row(k - 2);
	}

	//Lmd=(2*pi)*Lmd;
	//Lmd = Lmd * (2 * PI);
}

void UtilMath::polyleg(MatrixType& P, MatrixType x, unsigned n)
{
	//from polyleg.m

	if (x.cwiseAbs().maxCoeff() > 1) {
		cerr << "Values should be in range [-1, 1] \n";
		throw;
	}

	switch (n)
	{
	case 0:
		break;
	case 1:
		P.col(1) = x;
		break;
	default:
		P.col(1) = x;
		for (unsigned k = 2; k < n + 1; ++k)
		{
			double c1 = (2.0 * k - 1.0) / k;
			double c2 = (k - 1.0) / k;
			P.col(k) = c1 * x.cwiseProduct(P.col(k - 1)) - c2 * P.col(k - 2);
		}
	}
}
