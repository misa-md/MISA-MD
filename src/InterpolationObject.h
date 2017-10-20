class InterpolationObject{
public:
	InterpolationObject();
	~InterpolationObject();

	void initInterpolationObject(int _n, double _x0, double dx, double* data);
	void bcastInterpolationObject(int rank);
	void interpolatefile();
	int n;          //!< 表中数据个数
	double x0;      //!< 起始点
	double invDx;   //!< 倒数
	double* values;
	double** spline;
};