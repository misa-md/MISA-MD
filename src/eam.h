#ifndef EAM_H_
#define EAM_H_

#include <string>

#include "InterpolationObject.h"

/// 1 amu in kilograms
#define amuInKilograms  1.660538921e-27

/// 1 fs in seconds
#define fsInSeconds     1.0e-15

/// 1 Ang in meters
#define AngsInMeters    1.0e-10

/// 1 eV in Joules
#define eVInJoules      1.602176565e-19

/// Internal mass units are eV * fs^2 / Ang^2
static const double amuToInternalMass =
		amuInKilograms * AngsInMeters * AngsInMeters
		/ (fsInSeconds * fsInSeconds  * eVInJoules);

/// Hartrees to eVs
static const double hartreeToEv = 27.21138505;

/// Bohrs to Angstroms
static const double bohrToAngs = 0.52917721092;

using namespace std;

class eam{
public:
	eam();
	~eam();
	void setatomicNo(double nAtomic);
	void setlat(double latticeconst);
	void setmass(int i, double _mass);
	void setlatticeType(char* _latticeType);
	void setname(char* _name);
	void setcutoff(double _cutoff);
	void init(int nElems);
	void initf(int i, int nRho, double x0, double	dRho, double* buf);
	void initphi(int i, int nR, double x0, double	dR, double* buf);
	void initrho(int i, int nR, double x0, double	dR, double* buf);
	void eamBcast(int rank);

	void interpolatefile();
	InterpolationObject *phi;  //!< 对势
	InterpolationObject *rho;  //!< 电子云密度
	InterpolationObject *f;    //!< 嵌入能
private:
	int _nElems;
	double cutoff;          //!< 截断半径
	double *mass;           //!< 质量
	double lat;             //!< 晶格常数
	char latticeType[8];    //!< 晶格类型
	char  name[3];	   //!< 元素名
	int	 atomicNo;	   //!< 元素序号
};

#endif /* EAM_H_ */