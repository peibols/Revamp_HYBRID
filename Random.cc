#include <cmath>
#include "Random.h"

numrand::numrand()
{
}

numrand::numrand(int ir) : _ir(ir)
{
}

numrand::~numrand()
{
}

double numrand::rando()
{
        double da=16807.;
        double db=2147483647.;
        double dc=2147483648.;
        double usran;

	int use = _ir;
        _ir=int(fabs(fmod(da*use,db)+0.5));
        usran=double(_ir)/dc;
        return usran;
}

void numrand::SetIr(int ir)
{
	_ir=ir;
}
int numrand::GetIr() const
{
	return _ir;
}
