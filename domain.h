#ifndef DOMAIN_H_
#define DOMAIN_H_

class domain {
private:
	domain();
	domain(domain &domain);

	domain& operator=(domain &domain);
	
public:
	domain(int rank);
	double getGlobalLength(int index) const;
	void setGlobalLength(int index, double length);
	double getGlobalLength(int d) { return _globalLength[d]; }
private:
	int _localRank;
    double _globalLength[3];
};

#endif /*DOMAIN_H_*/
