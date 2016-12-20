#include "atom.h"
#include <iostream>
#include <math.h>
#include <fstream>
#include <iomanip>
#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836

atom::atom(double boxlo[3], double boxhi[3], double globalLengh[3], 
		   double boundingBoxMin[3], double boundingBoxMax[3], double ghostlength, double latticeconst, double cutoffRadius, int seed){
	_latticeconst = latticeconst;
	for (int d = 0; d < 3; d++) {
		_boxlo[d] = boxlo[d];
		_boxhi[d] = boxhi[d];
		_globalLengh[d] = globalLengh[d];
		_boundingBoxMin[d] = boundingBoxMin[d];
		_boundingBoxMax[d] = boundingBoxMax[d];
		_ghostLength[d] = ghostlength;
		_ghostBoundingBoxMin[d] = _boundingBoxMin[d] - _ghostLength[d];
		_ghostBoundingBoxMax[d] = _boundingBoxMax[d] + _ghostLength[d];
	}
	
	nlocalx = floor( _boundingBoxMax[0] / (_latticeconst) ) - floor( _boundingBoxMin[0] / (_latticeconst));
	nlocalx *= 2;
	nlocaly = floor( _boundingBoxMax[1] / _latticeconst ) - floor( _boundingBoxMin[1] / _latticeconst );
	nlocalz = floor( _boundingBoxMax[2] / _latticeconst ) - floor( _boundingBoxMin[2] / _latticeconst );

	/*
	nghostx = nlocalx + 2 * 2 * ( ceil( cutoffRadius / _latticeconst ) + 1 );
    nghosty = nlocaly + 2 * ( ceil( cutoffRadius / _latticeconst ) + 1 );
    nghostz = nlocalz + 2 * ( ceil( cutoffRadius / _latticeconst ) + 1 );
	*/

	nghostx = nlocalx + 2 * 2 * ceil( cutoffRadius / _latticeconst );
    nghosty = nlocaly + 2 * ceil( cutoffRadius / _latticeconst );
    nghostz = nlocalz + 2 * ceil( cutoffRadius / _latticeconst );


	lolocalx = floor( _boundingBoxMin[0] / latticeconst ) * 2;
	lolocaly = floor( _boundingBoxMin[1] / latticeconst );
	lolocalz = floor( _boundingBoxMin[2] / latticeconst );
	
	/*
	loghostx = lolocalx - 2 * ( ceil( cutoffRadius / _latticeconst ) + 1 );
    loghosty = lolocaly - ( ceil( cutoffRadius / _latticeconst ) + 1 );
    loghostz = lolocalz - ( ceil( cutoffRadius / _latticeconst ) + 1 );
	*/

	loghostx = lolocalx - 2 * ceil( cutoffRadius / _latticeconst );
    loghosty = lolocaly - ceil( cutoffRadius / _latticeconst );
    loghostz = lolocalz - ceil( cutoffRadius / _latticeconst );


	_seed = seed;
	_cutoffRadius = cutoffRadius;
	_cutlattice = ceil( _cutoffRadius / _latticeconst );

	nlocalinter = 0;
	nghostinter = 0;

	numberoflattice = nghostx * nghosty * nghostz;
	printf("number:%d, %d\n", numberoflattice, nlocalx*nlocaly*nlocalz);
	id = new unsigned long[numberoflattice];
	type = new int[numberoflattice];
	x = new double[numberoflattice * 3];
	v = new double[numberoflattice * 3];
	f = new double[numberoflattice * 3];
	rho = new double[numberoflattice];
	df = new double[numberoflattice];

	calculateNeighbourIndices();
}

atom::~atom(){
	delete[] id;
	delete[] type;
	delete[] x;
	delete[] v;
	delete[] f;
	delete[] rho;
	delete[] df;
}

void atom::calculateNeighbourIndices()
{
    int i;
    double x,y,z;
    int mark = 0;
	double cut_times_latti = _cutoffRadius / _latticeconst;
    vector<long int>::iterator neighbourOffsetsIter;
	for (int zIndex = -_cutlattice; zIndex <= _cutlattice; zIndex++)
	{
		for (int yIndex = -_cutlattice; yIndex <= _cutlattice; yIndex++)
		{
			for (int xIndex = -_cutlattice * 2; xIndex <= _cutlattice * 2; xIndex++)
			{
				//体心
				z=(double)zIndex+(((double)(xIndex%2))/2);
				y=(double)yIndex+(((double)(xIndex%2))/2);
				x=(double)xIndex/2;
				double square=x*x+y*y+z*z;
				long int offset;
				double r=sqrt(square);
				if (r<(cut_times_latti+0.4))
				{
					offset = IndexOf3DIndex(xIndex, yIndex, zIndex);
					if(offset > 0)
					{
						NeighbourOffsets.push_back(offset);
					} 
				} 

				//晶格点
				z=(double)zIndex-(((double)(xIndex%2))/2);
				y=(double)yIndex-(((double)(xIndex%2))/2);
				x=(double)xIndex/2;
				square = x*x+y*y+z*z;
				r=sqrt(square);
				if (r<(cut_times_latti+0.4))
				{
					offset = IndexOf3DIndex(xIndex, yIndex, zIndex);
					if(offset > 0)
					{
						for(neighbourOffsetsIter = NeighbourOffsets.begin(); neighbourOffsetsIter != NeighbourOffsets.end(); neighbourOffsetsIter++)
						{
							if( *neighbourOffsetsIter == offset)
							{
								mark = 1;                                
							}
						}
						if(mark != 1)
						{
							NeighbourOffsets.push_back(offset);
						}
						mark = 0;
					}
				} 
			} 
		}
	}
}

long int atom::IndexOf3DIndex(long int xIndex, long int yIndex, long int zIndex) const {
	return (zIndex * nghosty + yIndex) * nghostx + xIndex;
}

void atom::addatom(unsigned long id,double rx, double ry, double rz, double vx, double vy, double vz){
	int i;
	if ( ( rx >= _boundingBoxMin[0]) && (rx < _boundingBoxMax[0]) &&
		( ry >= _boundingBoxMin[1]) && (ry < _boundingBoxMax[1]) &&
		( rz >= _boundingBoxMin[2]) && (rz < _boundingBoxMax[2]) )
	{ 
		int lattice[3];
		lattice[0] = rx * 2 / _latticeconst + 0.5;
		lattice[1] = ry * 2 / _latticeconst + 0.5;
		lattice[2] = rz * 2 / _latticeconst + 0.5;
		lattice[1] = lattice[1] / 2;
		lattice[2] = lattice[2] / 2;
		lattice[0] -= loghostx;
		lattice[1] -= loghosty;
		lattice[2] -= loghostz;
		i = ((nghosty) * lattice[2] + lattice[1]) * (nghostx) + lattice[0];
		this->id[i] = id;
		x[i * 3] = rx;
		x[i * 3 + 1] = ry;
		x[i * 3 + 2] = rz;
		v[i * 3] = vx;
		v[i * 3 + 1] = vy;
		v[i * 3 + 2] = vz;
	}
}

int atom::decide(){
	nghostinter = 0;
	int nflag = 0;
	int kk = 0;
	double dist;
	double xtemp, ytemp, ztemp;
	int xstart = lolocalx - loghostx;
    int ystart = lolocaly - loghosty;
	int zstart = lolocalz - loghostz;
	
	//对本地晶格点原子进行判断，看是否运动为间隙原子
	for(int k = zstart; k < nlocalz + zstart; k++){
		for(int j = ystart; j < nlocaly + ystart; j++){
			for(int i = xstart; i < nlocalx + xstart; i++){
				kk = IndexOf3DIndex( i, j, k) * 3;
				if(x[kk] != -100){
					xtemp = (i +loghostx) * 0.5 * _latticeconst;
					ytemp = (j + loghosty + (i % 2) * 0.5) * _latticeconst;
					ztemp = (k + loghostz + (i % 2) * 0.5) * _latticeconst;
					dist = (x[kk] - xtemp) * (x[kk] - xtemp);
					dist += (x[kk + 1] - ytemp) * (x[kk + 1] - ytemp);
					dist += (x[kk + 2] - ztemp) * (x[kk + 2] - ztemp);
					if (dist > ( pow(0.2 * _latticeconst, 2.0) )/**超过距离则判断为间隙原子*/){
						if(xinter.size() > nlocalinter){
							if(idinter.size() > nlocalinter)
								idinter[nlocalinter] = id[kk / 3];
							else
								idinter.push_back(id[kk / 3]);
							typeinter[nlocalinter] = type[kk / 3];
							xinter[nlocalinter][0] = x[kk];
							xinter[nlocalinter][1] = x[kk + 1];
							xinter[nlocalinter][2] = x[kk + 2];
							if(vinter.size() > nlocalinter){
								vinter[nlocalinter][0] = v[kk];
								vinter[nlocalinter][1] = v[kk + 1];
								vinter[nlocalinter][2] = v[kk + 2];
							}
							else{
								vinter.resize(nlocalinter + 1, vector<double>(3));
								vinter[nlocalinter][0] = v[kk];
								vinter[nlocalinter][1] = v[kk + 1];
								vinter[nlocalinter][2] = v[kk + 2];
							}
							nlocalinter++;
							finter.resize(nlocalinter, vector<double>(3));
							rhointer.resize(nlocalinter);
							dfinter.resize(nlocalinter);
						}
						else{
							idinter.push_back(id[kk / 3]);
							typeinter.push_back(type[kk / 3]);
							xinter.resize(nlocalinter + 1, vector<double>(3));
							xinter[nlocalinter][0] = x[kk];
							xinter[nlocalinter][1] = x[kk + 1];
							xinter[nlocalinter][2] = x[kk + 2];
							vinter.resize(nlocalinter + 1, vector<double>(3));
							vinter[nlocalinter][0] = v[kk];
							vinter[nlocalinter][1] = v[kk + 1];
							vinter[nlocalinter][2] = v[kk + 2];
							nlocalinter++;
							finter.resize(nlocalinter, vector<double>(3));
							rhointer.resize(nlocalinter);
							dfinter.resize(nlocalinter);
						}

						x[kk] = -100;
						x[kk + 1] = -100;
						x[kk + 2] = -100;
						v[kk] = 0;
						v[kk + 1] = 0;
						v[kk + 2] = 0;
						nflag = 1;
					}
				}
			}
		}
	}
	for(int i = 0; i < nlocalinter; i++){
		if(xinter[i][0] < _boxlo[0]){
			xinter[i][0] += _globalLengh[0];
		}
		else if(xinter[i][0] >= _boxhi[0]){
			xinter[i][0] -= _globalLengh[0];
		}
		if(xinter[i][1] < _boxlo[1]){
			xinter[i][1] += _globalLengh[1];
		}
		else if(xinter[i][1] >= _boxhi[1]){
			xinter[i][1] -= _globalLengh[1];
		}
		if(xinter[i][2] < _boxlo[2]){
			xinter[i][2] += _globalLengh[2];
		}
		else if(xinter[i][2] >= _boxhi[2]){
			xinter[i][2] -= _globalLengh[2];
		}
	}

	//判断，如果跑出晶格点的?佑峙芑鼐Ц竦悖蚍呕鼐Ц竦闶榇娲⑵湫畔?
	for(int i = 0; i < nlocalinter; i++){
		int j, k, l;
		xtemp = xinter[i][0];
		ytemp = xinter[i][1];
		ztemp = xinter[i][2];
		j = xtemp * 2 / _latticeconst + 0.5;
		k = ytemp * 2 / _latticeconst + 0.5;
		l = ztemp * 2 / _latticeconst + 0.5;
		k = k / 2;
		l = l / 2;
		j -= loghostx;
		k -= loghosty;
		l -= loghostz;
		//判断是否在所表示晶格范围内
		if(j <= (nlocalx + 2 *( ceil( _cutoffRadius / _latticeconst ) + 1 )) 
			&& k <= (nlocaly + ( ceil( _cutoffRadius / _latticeconst ) + 1 )) 
			&& l <= (nlocalz + ( ceil( _cutoffRadius / _latticeconst ) + 1 ))){
				j = IndexOf3DIndex( j, k, l) * 3;
				if( x[j] == -100){
					id[j / 3] = idinter[i];
					type[j /3] = typeinter[i];
					x[j] = xinter[i][0];
					x[j + 1] = xinter[i][1];
					x[j + 2] = xinter[i][2];
					v[j] = vinter[i][0];
					v[j + 1] = vinter[i][1];
					v[j + 2] = vinter[i][2];

					idinter[i] = idinter[nlocalinter-1];
					typeinter[i] = typeinter[nlocalinter-1];
					xinter[i][0] = xinter[nlocalinter-1][0];
					xinter[i][1] = xinter[nlocalinter-1][1];
					xinter[i][2] = xinter[nlocalinter-1][2];
					vinter[i][0] = vinter[nlocalinter-1][0];
					vinter[i][1] = vinter[nlocalinter-1][1];
					vinter[i][2] = vinter[nlocalinter-1][2];

					i--;
					nlocalinter--;
				}
		}
	}
	return nflag;
}

void atom::clear_force(){
	for(int i = 0; i < numberoflattice * 3; i++){
		f[i] = 0;
	}
	for(int i = 0; i< numberoflattice; i++){
		rho[i] = 0;
	}
	for(int i = 0; i < finter.size(); i++){
		finter[i][0] = 0;
		finter[i][1] = 0;
		finter[i][2] = 0;
		rhointer[i] = 0;
	}
}

void atom::computeEam(eam* pot, domaindecomposition* _domaindecomposition, double& comm){
	double starttime, stoptime;
	double xtemp, ytemp, ztemp;
	double delx, dely, delz;
	vector<long int>::iterator neighbourOffsetsIter;
	InterpolationObject* rho_spline = &pot->rho[0];
	InterpolationObject* f_spline = &pot->f[0];
	InterpolationObject* phi_spline = &pot->phi[0];
	int n;
	double dist2;
	double r;
	double rhoTmp, dRho, dEmbed, dfEmbed, phiTmp, dPhi;
	int nr, m;
	double p;
	double** spline;
	double fpair;
	double recip, phi, phip, psip, z2, z2p;
    int kk;
	int xstart = lolocalx - loghostx;
	int ystart = lolocaly - loghosty;
	int zstart = lolocalz - loghostz;
	
	//本地晶格点上的原子计算电子云密度
	for(int k = zstart; k < nlocalz + zstart; k++){
		for(int j = ystart; j < nlocaly + ystart; j++){
			for(int i = xstart; i < nlocalx + xstart; i++){
				kk = IndexOf3DIndex( i, j, k);
				xtemp = x[kk * 3];
				ytemp = x[kk * 3 + 1];
				ztemp = x[kk * 3 + 2];
				if(xtemp != -100){
					//对晶格点邻居原子遍历
					for(neighbourOffsetsIter = NeighbourOffsets.begin(); neighbourOffsetsIter != NeighbourOffsets.end(); neighbourOffsetsIter++){
						n = (kk + *neighbourOffsetsIter);
						delx = xtemp - x[n * 3];
						dely = ytemp - x[n * 3 + 1];
						delz = ztemp - x[n * 3 + 2];
						dist2 = delx*delx + dely*dely + delz*delz;
						if (dist2 < ( _cutoffRadius * _cutoffRadius )){
							r = sqrt(dist2);
							nr = rho_spline->n;
							p = r * rho_spline->invDx + 1.0;
							m = static_cast<int> (p);
							m = std::max(1, std::min(m, (nr-1)));
							p -= m;
							p = std::min(p, 1.0);
							spline = rho_spline->spline;
							rhoTmp = ((spline[m][3]*p + spline[m][4])*p + spline[m][5])*p + spline[m][6];

							rho[kk] += rhoTmp;
							rho[n] += rhoTmp;
						}
					}
				}
			}
		}
	}

	//间隙原子电子云密度
	for(int i = 0; i < nlocalinter; i++){
		int j, k, l;
		xtemp = xinter[i][0];
		ytemp = xinter[i][1];
		ztemp = xinter[i][2];
		j = xtemp * 2 / _latticeconst + 0.5;
		k = ytemp * 2 / _latticeconst + 0.5;
		l = ztemp * 2 / _latticeconst + 0.5;
		k = k / 2;
		l = l / 2;
		j -= loghostx;
		k -= loghosty;
		l -= loghostz;
		j = IndexOf3DIndex( j, k, l);

		delx = xtemp - x[j * 3];
		dely = ytemp - x[j * 3 + 1];
		delz = ztemp - x[j * 3 + 2];
		dist2 = delx*delx + dely*dely + delz*delz;
		if (dist2 < ( _cutoffRadius * _cutoffRadius )){
			r = sqrt(dist2);
			nr = rho_spline->n;
			p = r * rho_spline->invDx + 1.0;
			m = static_cast<int> (p);
			m = std::max(1, std::min(m, (nr-1)));
			p -= m;
			p = std::min(p, 1.0);
			spline = rho_spline->spline;
			rhoTmp = ((spline[m][3]*p + spline[m][4])*p + spline[m][5])*p + spline[m][6];

			rhointer[i] += rhoTmp;
			rho[j] += rhoTmp;
		}
		for(neighbourOffsetsIter = NeighbourOffsets.begin(); neighbourOffsetsIter != NeighbourOffsets.end(); neighbourOffsetsIter++){
			//计算间隙原子的所有邻居
			n = (j + *neighbourOffsetsIter);
			delx = xtemp - x[n * 3];
			dely = ytemp - x[n * 3 + 1];
			delz = ztemp - x[n * 3 + 2];
			dist2 = delx*delx + dely*dely + delz*delz;
			if (dist2 < ( _cutoffRadius * _cutoffRadius )){
				r = sqrt(dist2);
				nr = rho_spline->n;
				p = r * rho_spline->invDx + 1.0;
				m = static_cast<int> (p);
				m = std::max(1, std::min(m, (nr-1)));
				p -= m;
				p = std::min(p, 1.0);
				spline = rho_spline->spline;
				rhoTmp = ((spline[m][3]*p + spline[m][4])*p + spline[m][5])*p + spline[m][6];

				rhointer[i] += rhoTmp;
				rho[n] += rhoTmp;
			}
			n = (j - *neighbourOffsetsIter);
			delx = xtemp - x[n * 3];
			dely = ytemp - x[n * 3 + 1];
			delz = ztemp - x[n * 3 + 2];
			dist2 = delx*delx + dely*dely + delz*delz;
			if (dist2 < ( _cutoffRadius * _cutoffRadius )){
				r = sqrt(dist2);
				nr = rho_spline->n;
				p = r * rho_spline->invDx + 1.0;
				m = static_cast<int> (p);
				m = std::max(1, std::min(m, (nr-1)));
				p -= m;
				p = std::min(p, 1.0);
				spline = rho_spline->spline;
				rhoTmp = ((spline[m][3]*p + spline[m][4])*p + spline[m][5])*p + spline[m][6];

				rhointer[i] += rhoTmp;
				rho[n] += rhoTmp;
			}
		}
		//对间隙原子遍历
		for(int k = i + 1; k < ( nghostinter + nlocalinter ); k++){
			delx = xtemp - xinter[k][0];
			dely = ytemp - xinter[k][1];
			delz = ztemp - xinter[k][2];
			dist2 = delx*delx + dely*dely + delz*delz;
			if (dist2 < ( _cutoffRadius * _cutoffRadius )){
				r = sqrt(dist2);
				nr = rho_spline->n;
				p = r * rho_spline->invDx + 1.0;
				m = static_cast<int> (p);
				m = std::max(1, std::min(m, (nr-1)));
				p -= m;
				p = std::min(p, 1.0);
				spline = rho_spline->spline;
				rhoTmp = ((spline[m][3]*p + spline[m][4])*p + spline[m][5])*p + spline[m][6];

				rhointer[i] += rhoTmp;
				rhointer[k] += rhoTmp;
			}
		}
		//计算间隙原子嵌入能导数
		nr = f_spline->n;
		p = rhointer[i] * f_spline->invDx + 1.0;
		m = static_cast<int> (p);
		m = std::max(1, std::min(m, (nr-1)));
		p -= m;
		p = std::min(p, 1.0);
		spline = f_spline->spline;
		dfEmbed = (spline[m][0]*p + spline[m][1])*p + spline[m][2];
		dfinter[i] = dfEmbed;
	}

	char tmp[20];
        sprintf(tmp, "rho.atom");
        ofstream outfile;
        /*outfile.open(tmp);
        for(int k = zstart; k < nlocalz + zstart; k++){
                for(int j = ystart; j < nlocaly + ystart; j++){
                        for(int i = xstart; i < nlocalx + xstart; i++){
                                kk = IndexOf3DIndex( i, j, k);
                                if(x[kk * 3] != -100)
                                        outfile << rho[kk] << std::endl;
                        }
                }
        }
	for(int i = 0; i < rho_spline->n; i++){
		outfile << rho_spline->spline[i][6] << std::endl;
	}
        outfile.close();*/

	//发送电子云密度
	starttime = MPI_Wtime();
	_domaindecomposition->sendrho(this);
	stoptime = MPI_Wtime();
        comm = stoptime - starttime;

        /*sprintf(tmp, "rho2.atom");
        outfile;
        outfile.open(tmp);
        for(int k = zstart; k < nlocalz + zstart; k++){
                for(int j = ystart; j < nlocaly + ystart; j++){
                        for(int i = xstart; i < nlocalx + xstart; i++){
                                kk = IndexOf3DIndex( i, j, k);
                                if(x[kk * 3] != -100)
                                        outfile << rho[kk] << std::endl;
                        }
                }
        }
        outfile.close();*/

	//本地晶格点计算嵌入能导数
	for(int k = zstart; k < nlocalz + zstart; k++){
		for(int j = ystart; j < nlocaly + ystart; j++){
			for(int i = xstart; i < nlocalx + xstart; i++){
				kk = IndexOf3DIndex( i, j, k);
				nr = f_spline->n;
				p = rho[kk] * f_spline->invDx + 1.0;
				m = static_cast<int> (p);
				m = std::max(1, std::min(m, (nr-1)));
				p -= m;
				p = std::min(p, 1.0);
				spline = f_spline->spline;
				dfEmbed = (spline[m][0]*p + spline[m][1])*p + spline[m][2];
				df[kk] = dfEmbed;
			}
		}
	}
	
	/*sprintf(tmp, "df.atom");
	outfile.open(tmp);
	for(int i = 0; i < f_spline->n; i++){
		outfile << i << " " << f_spline->spline[i][6] << std::endl;
	}
	outfile.close();*/

	//发送嵌入能导数
	starttime = MPI_Wtime();
	_domaindecomposition->sendDfEmbed(this);
	stoptime = MPI_Wtime();
        comm += stoptime - starttime;

	for(int k = zstart; k < nlocalz + zstart; k++){
		for(int j = ystart; j < nlocaly + ystart; j++){
			for(int i = xstart; i < nlocalx + xstart; i++){
				kk = IndexOf3DIndex( i, j, k);
				xtemp = x[kk * 3];
				ytemp = x[kk * 3 + 1];
				ztemp = x[kk * 3 + 2];
				if(xtemp != -100){
					//对晶格点邻居原子遍历
					for(neighbourOffsetsIter = NeighbourOffsets.begin(); neighbourOffsetsIter != NeighbourOffsets.end(); neighbourOffsetsIter++){
						n = (kk + *neighbourOffsetsIter);
						delx = xtemp - x[n * 3];
						dely = ytemp - x[n * 3 + 1];
						delz = ztemp - x[n * 3 + 2];
						dist2 = delx*delx + dely*dely + delz*delz;
						if (dist2 < ( _cutoffRadius * _cutoffRadius )){
							r = sqrt(dist2);
							nr = phi_spline->n;
							p = r * phi_spline->invDx + 1.0;
							m = static_cast<int> (p);
							m = std::max(1, std::min(m, (nr-1)));
							p -= m;
							p = std::min(p, 1.0);
							spline = phi_spline->spline;
							phiTmp = ((spline[m][3]*p + spline[m][4])*p + spline[m][5])*p + spline[m][6];
							dPhi = (spline[m][0]*p + spline[m][1])*p + spline[m][2];
							spline = rho_spline->spline;
							dRho = (spline[m][0]*p + spline[m][1])*p + spline[m][2];

							z2 = phiTmp;
							z2p = dPhi;
							recip = 1.0 / r;
							phi = z2 * recip;
							phip = z2p * recip - phi * recip;
							psip = (df[kk] + df[n])*dRho + phip;
							fpair = -psip * recip;

							f[kk * 3] += delx*fpair;
							f[kk * 3 + 1] += dely*fpair;
							f[kk * 3 + 2] += delz*fpair;

							f[n * 3] -= delx*fpair;
							f[n * 3 + 1] -= dely*fpair;
							f[n * 3 + 2] -= delz*fpair;
						}
					}
				}
			}
		}
	}

	/*sprintf(tmp, "f.atom");
        outfile.open(tmp);
        for(int i = 0; i < phi_spline->n; i++){
                outfile << i << " " << phi_spline->spline[i][6] << std::endl;
        }
        outfile.close();*/

	//间隙原子计算嵌入能和对势带来的力
	for(int i = 0; i < nlocalinter; i++){
		int j, k, l;
		xtemp = xinter[i][0];
		ytemp = xinter[i][1];
		ztemp = xinter[i][2];
		j = xtemp * 2 / _latticeconst + 0.5;
		k = ytemp * 2 / _latticeconst + 0.5;
		l = ztemp * 2 / _latticeconst + 0.5;
		k = k / 2;
		l = l / 2;
		j -= loghostx;
		k -= loghosty;
		l -= loghostz;
		j = IndexOf3DIndex( j, k, l);

		delx = xtemp - x[j * 3];
		dely = ytemp - x[j * 3 + 1];
		delz = ztemp - x[j * 3 + 2];
		dist2 = delx*delx + dely*dely + delz*delz;
		if (dist2 < ( _cutoffRadius * _cutoffRadius )){
			r = sqrt(dist2);
			nr = phi_spline->n;
			p = r * phi_spline->invDx + 1.0;
			m = static_cast<int> (p);
			m = std::max(1, std::min(m, (nr-1)));
			p -= m;
			p = std::min(p, 1.0);
			spline = phi_spline->spline;
			phiTmp = ((spline[m][3]*p + spline[m][4])*p + spline[m][5])*p + spline[m][6];
			dPhi = (spline[m][0]*p + spline[m][1])*p + spline[m][2];
			spline = rho_spline->spline;
			dRho = (spline[m][0]*p + spline[m][1])*p + spline[m][2];

			z2 = phiTmp;
			z2p = dPhi;
			recip = 1.0 / r;
			phi = z2 * recip;
			phip = z2p * recip - phi * recip;
			psip = (dfinter[i] + df[j])*dRho + phip;
			fpair = -psip * recip;

			finter[i][0] += delx*fpair;
			finter[i][1] += dely*fpair;
			finter[i][2] += delz*fpair;

			f[j * 3] -= delx*fpair;
			f[j * 3 + 1] -= dely*fpair;
			f[j * 3 + 2] -= delz*fpair;
		}
		for(neighbourOffsetsIter = NeighbourOffsets.begin(); neighbourOffsetsIter != NeighbourOffsets.end(); neighbourOffsetsIter++){
			//计算间隙原子的所有邻居
			n = (j + *neighbourOffsetsIter);
			delx = xtemp - x[n * 3];
			dely = ytemp - x[n * 3 + 1];
			delz = ztemp - x[n * 3 + 2];
			dist2 = delx*delx + dely*dely + delz*delz;
			if (dist2 < ( _cutoffRadius * _cutoffRadius )){
				r = sqrt(dist2);
				nr = phi_spline->n;
				p = r * phi_spline->invDx + 1.0;
				m = static_cast<int> (p);
				m = std::max(1, std::min(m, (nr-1)));
				p -= m;
				p = std::min(p, 1.0);
				spline = phi_spline->spline;
				phiTmp = ((spline[m][3]*p + spline[m][4])*p + spline[m][5])*p + spline[m][6];
				dPhi = (spline[m][0]*p + spline[m][1])*p + spline[m][2];
				spline = rho_spline->spline;
				dRho = (spline[m][0]*p + spline[m][1])*p + spline[m][2];

				z2 = phiTmp;
				z2p = dPhi;
				recip = 1.0 / r;
				phi = z2 * recip;
				phip = z2p * recip - phi * recip;
				psip = (dfinter[i] + df[n])*dRho + phip;
				fpair = -psip * recip;

				finter[i][0] += delx*fpair;
				finter[i][1] += dely*fpair;
				finter[i][2] += delz*fpair;

				f[n * 3] -= delx*fpair;
				f[n * 3 + 1] -= dely*fpair;
				f[n * 3 + 2] -= delz*fpair;
			}
			n = (j - *neighbourOffsetsIter);
			delx = xtemp - x[n * 3];
			dely = ytemp - x[n * 3 + 1];
			delz = ztemp - x[n * 3 + 2];
			dist2 = delx*delx + dely*dely + delz*delz;
			if (dist2 < ( _cutoffRadius * _cutoffRadius )){
				r = sqrt(dist2);
				nr = phi_spline->n;
				p = r * phi_spline->invDx + 1.0;
				m = static_cast<int> (p);
				m = std::max(1, std::min(m, (nr-1)));
				p -= m;
				p = std::min(p, 1.0);
				spline = phi_spline->spline;
				phiTmp = ((spline[m][3]*p + spline[m][4])*p + spline[m][5])*p + spline[m][6];
				dPhi = (spline[m][0]*p + spline[m][1])*p + spline[m][2];
				spline = rho_spline->spline;
				dRho = (spline[m][0]*p + spline[m][1])*p + spline[m][2];

				z2 = phiTmp;
				z2p = dPhi;
				recip = 1.0 / r;
				phi = z2 * recip;
				phip = z2p * recip - phi * recip;
				psip = (dfinter[i] + df[n])*dRho + phip;
				fpair = -psip * recip;

				finter[i][0] += delx*fpair;
				finter[i][1] += dely*fpair;
				finter[i][2] += delz*fpair;

				f[n * 3] -= delx*fpair;
				f[n * 3 + 1] -= dely*fpair;
				f[n * 3 + 2] -= delz*fpair;
			}
		}
		//对间隙原子遍历
		for(int k = i + 1; k < ( nghostinter + nlocalinter ); k++){
			delx = xtemp - xinter[k][0];
			dely = ytemp - xinter[k][1];
			delz = ztemp - xinter[k][2];
			dist2 = delx*delx + dely*dely + delz*delz;
			if (dist2 < ( _cutoffRadius * _cutoffRadius )){
				r = sqrt(dist2);
				nr = phi_spline->n;
				p = r * phi_spline->invDx + 1.0;
				m = static_cast<int> (p);
				m = std::max(1, std::min(m, (nr-1)));
				p -= m;
				p = std::min(p, 1.0);
				spline = phi_spline->spline;
				phiTmp = ((spline[m][3]*p + spline[m][4])*p + spline[m][5])*p + spline[m][6];
				dPhi = (spline[m][0]*p + spline[m][1])*p + spline[m][2];
				spline = rho_spline->spline;
				dRho = (spline[m][0]*p + spline[m][1])*p + spline[m][2];

				z2 = phiTmp;
				z2p = dPhi;
				recip = 1.0 / r;
				phi = z2 * recip;
				phip = z2p * recip - phi * recip;
				psip = (dfinter[i] + dfinter[k])*dRho + phip;
				fpair = -psip * recip;

				finter[i][0] += delx*fpair;
				finter[i][1] += dely*fpair;
				finter[i][2] += delz*fpair;

				finter[k][0] -= delx*fpair;
				finter[k][1] -= dely*fpair;
				finter[k][2] -= delz*fpair;
			}
		}
	}
}

double atom::getBoundingBoxMin(int dimension) const {
	return this->_boundingBoxMin[dimension];
}

double atom::getBoundingBoxMax(int dimension) const {
	return this->_boundingBoxMax[dimension];
}

double atom::get_ghostlengh(int index) const {
	return _ghostLength[index];
}

int atom::getinteridsendsize(){
	return interbuf.size();
}

void atom::getatomx(int direction, vector<vector<int> > &sendlist){
	int i;
	if(direction == 0){
		//找到要发送到邻居进程的区域
		int xstart = lolocalx - loghostx;
		int ystart = lolocaly - loghosty;
		int zstart = lolocalz - loghostz;
		int xstop = xstart + (_cutlattice) * 2;
		int ystop = ystart + nlocaly;
		int zstop = zstart + nlocalz;
		
		//要发送要邻居进程区域内的分子指针
		for (int iz = zstart; iz < zstop; iz++) {
			for (int iy = ystart; iy < ystop; iy++) {
				for (int ix = xstart; ix < xstop; ix++) {
					i = IndexOf3DIndex( ix, iy, iz);
					sendlist[0].push_back(i);
				}
			}
		}
	}
	else{
		//找到要发送到邻居进程的区域
		int xstart = lolocalx - loghostx + nlocalx - ((_cutlattice) * 2);
		int ystart = lolocaly - loghosty;
		int zstart = lolocalz - loghostz;
		int xstop = lolocalx - loghostx + nlocalx;
		int ystop = ystart + nlocaly;
		int zstop = zstart + nlocalz;
        
		//要发送要邻居进程区域内的分子指针
		for (int iz = zstart; iz < zstop; iz++) {
			for (int iy = ystart; iy < ystop; iy++) {
				for (int ix = xstart; ix < xstop; ix++) {
					i = IndexOf3DIndex( ix, iy, iz);
					sendlist[1].push_back(i);
				}
			}
		}
	}
}

void atom::getatomy(int direction, vector<vector<int> > &sendlist){
	int i;
	if(direction == 0){
		//找到要发送到邻居进程的区域
		int xstart = 0;
		int ystart = lolocaly - loghosty;
		int zstart = lolocalz - loghostz;
		int xstop = nghostx;
		int ystop = ystart + _cutlattice;
		int zstop = zstart + nlocalz;
		
		//要发送要邻居进程区域内的分子指针
		for (int iz = zstart; iz < zstop; iz++) {
			for (int iy = ystart; iy < ystop; iy++) {
				for (int ix = xstart; ix < xstop; ix++) {
					i = IndexOf3DIndex( ix, iy, iz);
					sendlist[2].push_back(i);
				}
			}
		}
	}
	else{
		//找到要发送到邻居进程的区域
		int xstart = 0;
		int ystart = lolocaly - loghosty + nlocaly - (_cutlattice);
		int zstart = lolocalz - loghostz;
		int xstop = nghostx;
		int ystop = lolocaly - loghosty + nlocaly;
		int zstop = zstart + nlocalz;
        
		//要发送要邻居进程区域内的分子指针
		for (int iz = zstart; iz < zstop; iz++) {
			for (int iy = ystart; iy < ystop; iy++) {
				for (int ix = xstart; ix < xstop; ix++) {
					i = IndexOf3DIndex( ix, iy, iz);
					sendlist[3].push_back(i);
				}
			}
		}
	}
}

void atom::getatomz(int direction, vector<vector<int> > &sendlist){
	int i;
	if(direction == 0){
		//找到要发送到邻居进程的区域
		int xstart = 0;
		int ystart = 0;
		int zstart = lolocalz - loghostz;
		int xstop = nghostx;
		int ystop = nghosty;
		int zstop = zstart + _cutlattice;

		//要发送要邻居进程区域内的分子指针
		for (int iz = zstart; iz < zstop; iz++) {
			for (int iy = ystart; iy < ystop; iy++) {
				for (int ix = xstart; ix < xstop; ix++) {
					i = IndexOf3DIndex( ix, iy, iz);
					sendlist[4].push_back(i);
				}
			}
		}
	}
	else{
		//找到要发送到邻居进程的区域
		int xstart = 0;
		int ystart = 0;
		int zstart = lolocalz - loghostz + nlocalz - (_cutlattice);
		int xstop = nghostx;
		int ystop = nghosty;
		int zstop = lolocalz - loghostz + nlocalz;

		//要发送要邻居进程区域内的分子指针
		for (int iz = zstart; iz < zstop; iz++) {
			for (int iy = ystart; iy < ystop; iy++) {
				for (int ix = xstart; ix < xstop; ix++) {
					i = IndexOf3DIndex( ix, iy, iz);
					sendlist[5].push_back(i);
				}
			}
		}
	}
}

void atom::getIntertosend(int d, int direction, double ghostlengh, vector<int> &sendlist){
	double low, high;
	if(d == 0){
		if(direction == 0){
			low = _boundingBoxMin[0];
			high = _boundingBoxMin[0] + ghostlengh;
			for(int i = 0; i < nlocalinter; i++){
				if(xinter[i][0] < high && xinter[i][0] >= low){
					sendlist.push_back(i);
				}
			}
		}
		else{
			low = _boundingBoxMax[0] - ghostlengh;
			high = _boundingBoxMax[0];
			for(int i = 0; i < nlocalinter; i++){
				if(xinter[i][0] <= high && xinter[i][0] > low){
					sendlist.push_back(i);
				}
			}
		}
	}
	else if(d == 1){
		if(direction == 0){
			low = _boundingBoxMin[1];
			high = _boundingBoxMin[1] + ghostlengh;
			for(int i = 0; i < nlocalinter + nghostinter; i++){
				if(xinter[i][1] < high && xinter[i][1] >= low){
					sendlist.push_back(i);
				}
			}
		}
		else{
			low = _boundingBoxMax[1] - ghostlengh;
			high = _boundingBoxMax[1];
			for(int i = 0; i < nlocalinter + nghostinter; i++){
				if(xinter[i][1] <= high && xinter[i][1] > low){
					sendlist.push_back(i);
				}
			}
		}
	}
	else{
		if(direction == 0){
			low = _boundingBoxMin[2];
			high = _boundingBoxMin[2] + ghostlengh;
			for(int i = 0; i < nlocalinter + nghostinter; i++){
				if(xinter[i][2] < high && xinter[i][2] >= low){
					sendlist.push_back(i);
				}
			}
		}
		else{
			low = _boundingBoxMax[2] - ghostlengh;
			high = _boundingBoxMax[2];
			for(int i = 0; i < nlocalinter + nghostinter; i++){
				if(xinter[i][2] <= high && xinter[i][2] > low){
					sendlist.push_back(i);
				}
			}
		}
	}
}

int atom::getintersendnum(int dimension, int direction){
	interbuf.clear();
	for(int i = 0; i < nlocalinter; i++){
		if(direction == 0){
			if(xinter[i][dimension] < _boundingBoxMin[dimension]){
				interbuf.push_back(i);
			}
		}
		else{
			if(xinter[i][dimension] >= _boundingBoxMax[dimension]){
				interbuf.push_back(i);
			}
		}
	}
	return interbuf.size();
}

void atom::pack_intersend(particledata *buf){
	int j;
	for(int i = 0; i < interbuf.size(); i++){
		    j = interbuf[i];
			buf[i].id = idinter[j];
			buf[i].type = typeinter[j];
			buf[i].r[0] = xinter[j][0];
			buf[i].r[1] = xinter[j][1];
			buf[i].r[2] = xinter[j][2];
			buf[i].v[0] = vinter[j][0];
			buf[i].v[1] = vinter[j][1];
			buf[i].v[2] = vinter[j][2];
			idinter[j] = idinter[nlocalinter-1];
			typeinter[j] = typeinter[nlocalinter-1];
			xinter[j][0] = xinter[nlocalinter-1][0];
			xinter[j][1] = xinter[nlocalinter-1][1];
			xinter[j][2] = xinter[nlocalinter-1][2];
			vinter[j][0] = vinter[nlocalinter-1][0];
			vinter[j][1] = vinter[nlocalinter-1][1];
			vinter[j][2] = vinter[nlocalinter-1][2];
			nlocalinter--;
	}
}

void atom::unpack_interrecv(int d, int n, particledata *buf){
	vector<double> xtemp(3);
	vector<double> vtemp(3);
	unsigned long id;
	int type;
	int m = 0;
	for(int i = 0; i < n; i++){
		id = buf[i].id;
		type = buf[i].type;
		xtemp[0] = buf[i].r[0];
		xtemp[1] = buf[i].r[1];
		xtemp[2] = buf[i].r[2];
		vtemp[0] = buf[i].v[0];
		vtemp[1] = buf[i].v[1];
		vtemp[2] = buf[i].v[2];
		if(xtemp[d] >= _boundingBoxMin[d] && xtemp[d] < _boundingBoxMax[d]){
			if(nlocalinter == xinter.size()){
				idinter.push_back(id);
				typeinter.push_back(type);
				xinter.push_back(xtemp);
				vinter.push_back(vtemp);
				nlocalinter++;
				finter.resize(nlocalinter, vector<double>(3));
				rhointer.resize(nlocalinter);
				dfinter.resize(nlocalinter);
			}
			else{
				if(idinter.size() == nlocalinter)
					idinter.push_back(id);
				else
					idinter[nlocalinter] = id;
				typeinter[nlocalinter] = type;
				xinter[nlocalinter][0] = xtemp[0];
				xinter[nlocalinter][1] = xtemp[1];
				xinter[nlocalinter][2] = xtemp[2];
				if(nlocalinter == vinter.size()){
					vinter.push_back(vtemp);
				}
				else{
					vinter[nlocalinter][0] = vtemp[0];
					vinter[nlocalinter][1] = vtemp[1];
					vinter[nlocalinter][2] = vtemp[2];
				}
				nlocalinter++;
				finter.resize(nlocalinter, vector<double>(3));
				rhointer.resize(nlocalinter);
				dfinter.resize(nlocalinter);
			}
		}
	}
}

void atom::pack_bordersend(int dimension, int n, vector<int> &sendlist, latparticledata *buf, double shift){
	int j;
	if(dimension == 0){
		for(int i = 0; i < n; i++){
			j = sendlist[i];
			buf[i].type = typeinter[j];
			buf[i].r[0] = xinter[j][0] + shift;
			buf[i].r[1] = xinter[j][1];
			buf[i].r[2] = xinter[j][2];
		}
	}
	else if(dimension == 1){
		for(int i = 0; i < n; i++){
			j = sendlist[i];
			buf[i].type = typeinter[j];
			buf[i].r[0] = xinter[j][0];
			buf[i].r[1] = xinter[j][1] + shift;
			buf[i].r[2] = xinter[j][2];
		}
	}
	else{
		for(int i = 0; i < n; i++){
			j = sendlist[i];
			buf[i].type = typeinter[j];
			buf[i].r[0] = xinter[j][0];
			buf[i].r[1] = xinter[j][1];
			buf[i].r[2] = xinter[j][2] + shift;
		}
	}
}

void atom::unpack_borderrecv(int n, latparticledata *buf, vector<int> &recvlist){
	int type;
	vector<double> xtemp(3);
	for(int i = 0; i < n; i++){
		type = buf[i].type;
		xtemp[0] = buf[i].r[0];
		xtemp[1] = buf[i].r[1];
		xtemp[2] = buf[i].r[2];
		if(xtemp[0] >= _ghostBoundingBoxMin[0] && xtemp[0] < _ghostBoundingBoxMax[0] &&
			xtemp[1] >= _ghostBoundingBoxMin[1] && xtemp[1] < _ghostBoundingBoxMax[1] &&
			xtemp[2] >= _ghostBoundingBoxMin[2] && xtemp[2] < _ghostBoundingBoxMax[2]){
				if(xinter.size() == nlocalinter + nghostinter){
					typeinter.push_back(type);
					xinter.push_back(xtemp);
					nghostinter++;
					recvlist[i] = nlocalinter + nghostinter - 1;
					finter.resize(nlocalinter + nghostinter, vector<double>(3));
					rhointer.resize(nlocalinter + nghostinter);
					dfinter.resize(nlocalinter + nghostinter);
				}
				else{
					typeinter[nlocalinter + nghostinter] = type;
					xinter[nlocalinter + nghostinter][0] = xtemp[0];
					xinter[nlocalinter + nghostinter][1] = xtemp[1];
					xinter[nlocalinter + nghostinter][2] = xtemp[2];
					nghostinter++;
					recvlist[i] = nlocalinter + nghostinter - 1;
					finter.resize(nlocalinter + nghostinter, vector<double>(3));
					rhointer.resize(nlocalinter + nghostinter);
					dfinter.resize(nlocalinter + nghostinter);
				}
		}
		else
			recvlist[i] = -1;
	}
}

void atom::pack_send(int dimension, int n, vector<int> &sendlist, latparticledata *buf, double shift){
	int j;
	if(dimension == 0){
		for(int i = 0; i < n; i++){
			j = sendlist[i];
			buf[i].type = type[j];
			buf[i].r[0] = x[j * 3] + shift;
			buf[i].r[1] = x[j * 3 + 1];
			buf[i].r[2] = x[j * 3 + 2];
		}
	}
	else if(dimension == 1){
		for(int i = 0; i < n; i++){
			j = sendlist[i];
			buf[i].type = type[j];
			buf[i].r[0] = x[j * 3];
			buf[i].r[1] = x[j * 3 + 1] + shift;
			buf[i].r[2] = x[j * 3 + 2];
		}
	}
	else{
		for(int i = 0; i < n; i++){
			j = sendlist[i];
			buf[i].type = type[j];
			buf[i].r[0] = x[j * 3];
			buf[i].r[1] = x[j * 3 + 1];
			buf[i].r[2] = x[j * 3 + 2] + shift;
		}
	}
}

void atom::unpack_recvfirst(int d, int direction, latparticledata *buf, vector<vector<int> > &recvlist){
	int xstart, ystart, zstart;
	int xstop, ystop, zstop;
	int kk;
	int m = 0;
	if(d == 0){
		if(direction == 0){
			xstart = lolocalx - loghostx + nlocalx;
			xstop = nghostx;
			ystart = lolocaly - loghosty;
			ystop = ystart + nlocaly;
			zstart = lolocalz - loghostz;
			zstop = zstart + nlocalz;
			for(int k = zstart; k < zstop; k++){
				for(int j = ystart; j < ystop; j++){
					for(int i = xstart; i < xstop; i++){
						kk = IndexOf3DIndex(i, j, k);
						type[kk] = buf[m].type;
						x[kk * 3] = buf[m].r[0];
						x[kk * 3 + 1] = buf[m].r[1];
						x[kk * 3 + 2] = buf[m++].r[2];
						recvlist[0].push_back(kk);
					}
				}
			}
		}
		else{
			xstart = 0;
			xstop = lolocalx - loghostx;
			ystart = lolocaly - loghosty;
			ystop = ystart + nlocaly;
			zstart = lolocalz - loghostz;
			zstop = zstart + nlocalz;
			for(int k = zstart; k < zstop; k++){
				for(int j = ystart; j < ystop; j++){
					for(int i = xstart; i < xstop; i++){
						kk = IndexOf3DIndex(i, j, k);
						type[kk] = buf[m].type;
						x[kk * 3] = buf[m].r[0];
						x[kk * 3 + 1] = buf[m].r[1];
						x[kk * 3 + 2] = buf[m++].r[2];
						recvlist[1].push_back(kk);
					}
				}
			}
		}
	}
	else if (d == 1){
		if(direction == 0){
			xstart = 0;
			xstop = nghostx;
			ystart = lolocaly - loghosty + nlocaly;
			ystop = nghosty;
			zstart = lolocalz - loghostz;
			zstop = zstart + nlocalz;
			for(int k = zstart; k < zstop; k++){
				for(int j = ystart; j < ystop; j++){
					for(int i = xstart; i < xstop; i++){
						kk = IndexOf3DIndex(i, j, k);
						type[kk] = buf[m].type;
						x[kk * 3] = buf[m].r[0];
						x[kk * 3 + 1] = buf[m].r[1];
						x[kk * 3 + 2] = buf[m++].r[2];
						recvlist[2].push_back(kk);
					}
				}
			}
		}
		else{
			xstart = 0;
			xstop = nghostx;
			ystart = 0;
			ystop = lolocaly - loghosty;
			zstart = lolocalz - loghostz;
			zstop = zstart + nlocalz;
			for(int k = zstart; k < zstop; k++){
				for(int j = ystart; j < ystop; j++){
					for(int i = xstart; i < xstop; i++){
						kk = IndexOf3DIndex(i, j, k);
						type[kk] = buf[m].type;
						x[kk * 3] = buf[m].r[0];
						x[kk * 3 + 1] = buf[m].r[1];
						x[kk * 3 + 2] = buf[m++].r[2];
						recvlist[3].push_back(kk);
					}
				}
			}
		}
	}
	else{
		if(direction == 0){
			xstart = 0;
			xstop = nghostx;
			ystart = 0;
			ystop = nghosty;
			zstart = lolocalz - loghostz + nlocalz;
			zstop = nghostz;
			for(int k = zstart; k < zstop; k++){
				for(int j = ystart; j < ystop; j++){
					for(int i = xstart; i < xstop; i++){
						kk = IndexOf3DIndex(i, j, k);
						type[kk] = buf[m].type;
						x[kk * 3] = buf[m].r[0];
						x[kk * 3 + 1] = buf[m].r[1];
						x[kk * 3 + 2] = buf[m++].r[2];
						recvlist[4].push_back(kk);
					}
				}
			}
		}
		else{
			xstart = 0;
			xstop = nghostx;
			ystart = 0;
			ystop = nghosty;
			zstart = 0;
			zstop = lolocalz - loghostz;
			for(int k = zstart; k < zstop; k++){
				for(int j = ystart; j < ystop; j++){
					for(int i = xstart; i < xstop; i++){
						kk = IndexOf3DIndex(i, j, k);
						type[kk] = buf[m].type;
						x[kk * 3] = buf[m].r[0];
						x[kk * 3 + 1] = buf[m].r[1];
						x[kk * 3 + 2] = buf[m++].r[2];
						recvlist[5].push_back(kk);
					}
				}
			}
		}
	}
}

void atom::unpack_recv(int d, int direction, int n, latparticledata *buf, vector<vector<int> > &recvlist){
	int kk;
	if(d == 0){
		if(direction == 0){
			for(int i = 0; i < n; i++){
				kk = recvlist[0][i];
				type[kk] = buf[i].type;
				x[kk * 3] = buf[i].r[0];
				x[kk * 3 + 1] = buf[i].r[1];
				x[kk * 3 + 2] = buf[i].r[2];
			}
		}
		else{
			for(int i = 0; i < n; i++){
				kk = recvlist[1][i];
				type[kk] = buf[i].type;
				x[kk * 3] = buf[i].r[0];
				x[kk * 3 + 1] = buf[i].r[1];
				x[kk * 3 + 2] = buf[i].r[2];
			}
		}
	}
	else if (d == 1){
		if(direction == 0){
			for(int i = 0; i < n; i++){
				kk = recvlist[2][i];
				type[kk] = buf[i].type;
				x[kk * 3] = buf[i].r[0];
				x[kk * 3 + 1] = buf[i].r[1];
				x[kk * 3 + 2] = buf[i].r[2];
			}
		}
		else{
			for(int i = 0; i < n; i++){
				kk = recvlist[3][i];
				type[kk] = buf[i].type;
				x[kk * 3] = buf[i].r[0];
				x[kk * 3 + 1] = buf[i].r[1];
				x[kk * 3 + 2] = buf[i].r[2];
			}
		}
	}
	else{
		if(direction == 0){
			for(int i = 0; i < n; i++){
				kk = recvlist[4][i];
				type[kk] = buf[i].type;
				x[kk * 3] = buf[i].r[0];
				x[kk * 3 + 1] = buf[i].r[1];
				x[kk * 3 + 2] = buf[i].r[2];
			}
		}
		else{
			for(int i = 0; i < n; i++){
				kk = recvlist[5][i];
				type[kk] = buf[i].type;
				x[kk * 3] = buf[i].r[0];
				x[kk * 3 + 1] = buf[i].r[1];
				x[kk * 3 + 2] = buf[i].r[2];
			}
		}
	}
}

void atom::pack_rho(int n, vector<int> &recvlist, double *buf){
	int j, m = 0;
	for(int i = 0; i < n; i++){
		j = recvlist[i];
		buf[m++] = rho[j];
	}
}

void atom::unpack_rho(int d, int direction, double *buf, vector<vector<int> > &sendlist){
	int j, m = 0;
	if(d == 0){
		if(direction == 0){
			for(int i = 0; i < sendlist[1].size(); i++){
				j = sendlist[1][i];
				rho[j] += buf[m++];
			}
		}
		else{
			for(int i = 0; i < sendlist[0].size(); i++){
				j = sendlist[0][i];
				rho[j] += buf[m++];
			}
		}
	}
	else if(d == 1){
		if(direction == 0){
			for(int i = 0; i < sendlist[3].size(); i++){
				j = sendlist[3][i];
				rho[j] += buf[m++];
			}
		}
		else{
			for(int i = 0; i < sendlist[2].size(); i++){
				j = sendlist[2][i];
				rho[j] += buf[m++];
			}
		}
	}
	else{
		if(direction == 0){
			for(int i = 0; i < sendlist[5].size(); i++){
				j = sendlist[5][i];
				rho[j] += buf[m++];
			}
		}
		else{
			for(int i = 0; i < sendlist[4].size(); i++){
				j = sendlist[4][i];
				rho[j] += buf[m++];
			}
		}
	}
}

void atom::pack_df(vector<int> &sendlist, vector<int> &intersendlist, double *buf){
	int j, m = 0;
	int n = sendlist.size();
	for(int i = 0; i < n; i++){
		j = sendlist[i];
		buf[m++] = df[j];
	}
	n = intersendlist.size();
	for(int i = 0; i < n; i++){
		j = intersendlist[i];
		buf[m++] = dfinter[j];
	}
}

void atom::unpack_df(int n, double *buf, vector<int> &recvlist, vector<int> &interrecvlist){
	int kk;
	int m = 0;
	if(n != ( recvlist.size() + interrecvlist.size() )){
		printf("wrong number of dfembed recv!!!");
		MPI_Abort(MPI_COMM_WORLD, 2);
	}
	n = recvlist.size();
	for(int i = 0; i < n; i++){
		kk = recvlist[i];
		df[kk] = buf[m++];
	}
	n = interrecvlist.size();
	for(int i = 0; i < n; i++){
		kk = interrecvlist[i];
		if(kk == -1){
			m++;
			continue;
		}
		dfinter[kk] = buf[m++];
	}
}

void atom::pack_force(int n, vector<int> &recvlist, double *buf){
	int j, m = 0;
	for(int i = 0; i < n; i++){
		j = recvlist[i] * 3;
		buf[m++] = f[j];
		buf[m++] = f[j + 1];
		buf[m++] = f[j + 2];
	}
}

void atom::unpack_force(int d, int direction, double *buf, vector<vector<int> > &sendlist){
	int j, m = 0;
	if(d == 0){
		if(direction == 0){
			for(int i = 0; i < sendlist[1].size(); i++){
				j = sendlist[1][i] * 3;
				f[j] += buf[m++];
				f[j + 1] += buf[m++];
				f[j + 2] += buf[m++];
			}
		}
		else{
			for(int i = 0; i < sendlist[0].size(); i++){
				j = sendlist[0][i] * 3;
				f[j] += buf[m++];
				f[j + 1] += buf[m++];
				f[j + 2] += buf[m++];
			}
		}
	}
	else if(d == 1){
		if(direction == 0){
			for(int i = 0; i < sendlist[3].size(); i++){
				j = sendlist[3][i] * 3;
				f[j] += buf[m++];
				f[j + 1] += buf[m++];
				f[j + 2] += buf[m++];
			}
		}
		else{
			for(int i = 0; i < sendlist[2].size(); i++){
				j = sendlist[2][i] * 3;
				f[j] += buf[m++];
				f[j + 1] += buf[m++];
				f[j + 2] += buf[m++];
			}
		}
	}
	else{
		if(direction == 0){
			for(int i = 0; i < sendlist[5].size(); i++){
				j = sendlist[5][i] * 3;
				f[j] += buf[m++];
				f[j + 1] += buf[m++];
				f[j + 2] += buf[m++];
			}
		}
		else{
			for(int i = 0; i < sendlist[4].size(); i++){
				j = sendlist[4][i] * 3;
				f[j] += buf[m++];
				f[j + 1] += buf[m++];
				f[j + 2] += buf[m++];
			}
		}
	}
}

void atom::computefirst(double dtInv2m, double dt){
	int kk;
	int xstart = lolocalx - loghostx;
	int ystart = lolocaly - loghosty;
	int zstart = lolocalz - loghostz;
	//本地晶格点上的原子求解运动方程第一步
	for(int k = zstart; k < nlocalz + zstart; k++){
		for(int j = ystart; j < nlocaly + ystart; j++){
			for(int i = xstart; i < nlocalx + xstart; i++){
				kk = IndexOf3DIndex( i, j, k) * 3;
				if(x[kk] != -100){
					v[kk] = v[kk] + dtInv2m * f[kk];
					x[kk] += dt * v[kk];
					v[kk + 1] = v[kk + 1] + dtInv2m * f[kk + 1];
					x[kk + 1] += dt * v[kk + 1];
					v[kk + 2] = v[kk + 2] + dtInv2m * f[kk + 2];
					x[kk + 2] += dt * v[kk + 2];
				}
			}
		}
	}
    //本地间隙原子求解运动方程第一步
	for(int i = 0; i < nlocalinter; i++){
		for (unsigned short d = 0; d < 3; ++d) {
			vinter[i][d] = vinter[i][d] + dtInv2m * finter[i][d];
			xinter[i][d] += dt * vinter[i][d];
		}
	}
}

void atom::computesecond(double dtInv2m){
	int kk;
	int xstart = lolocalx - loghostx;
	int ystart = lolocaly - loghosty;
	int zstart = lolocalz - loghostz;
	//本地晶格点上的原子求解运动方程第二步
	for(int k = zstart; k < nlocalz + zstart; k++){
		for(int j = ystart; j < nlocaly + ystart; j++){
			for(int i = xstart; i < nlocalx + xstart; i++){
				kk = IndexOf3DIndex( i, j, k) * 3;
				for (unsigned short d = 0; d < 3; ++d) {
					v[kk + d] += dtInv2m * f[kk + d];
				}
			}
		}
	}
    //本地间隙原子求解运动方程第二步
	for(int i = 0; i < nlocalinter; i++){
		for (unsigned short d = 0; d < 3; ++d) {
			vinter[i][d] = vinter[i][d] + dtInv2m * finter[i][d];
		}
	}
}

void atom::print_force(){
        char tmp[20];
        sprintf(tmp, "force.txt");
        ofstream outfile;
        outfile.open(tmp);
        int kk;
        int xstart = lolocalx - loghostx;
        int ystart = lolocaly - loghosty;
        int zstart = lolocalz - loghostz;
        std::cout << "print_force" << std::endl;
        for(int k = zstart; k < nlocalz + zstart; k++){
                for(int j = ystart; j < nlocaly + ystart; j++){
                        for(int i = xstart; i < nlocalx + xstart; i++){
                                kk = IndexOf3DIndex( i, j, k) * 3;
                                outfile << f[kk] << " " << f[kk + 1] << " " << f[kk + 2] << std::endl;
                        }
                }
        }
        outfile.close();
}

void atom::setv(int lat[4], double collision_v[3]){
	int kk;
	if((lat[0] * 2) >= lolocalx && (lat[0] * 2) < (lolocalx + nlocalx)
		&& lat[1] >= lolocaly && lat[1] < (lolocaly + nlocaly)
		&& lat[2] >= lolocalz && lat[2] < (lolocalz + nlocalz))
	{
		kk = (IndexOf3DIndex(lat[0] * 2-loghostx, lat[1]-loghosty, lat[2]-loghostz) + lat[4]) * 3;
		v[kk] += collision_v[0];
		v[kk + 1] += collision_v[1];
		v[kk + 2] += collision_v[2];
	}
}

void atom::print_atom(int rank){
	char tmp[20];
	sprintf(tmp, "dump_%d.atom", rank);
	ofstream outfile;
	outfile.open(tmp);
	int kk;
	int xstart = lolocalx - loghostx;
	int ystart = lolocaly - loghosty;
	int zstart = lolocalz - loghostz;
	outfile << "print_atom" << std::endl;
	for(int k = zstart; k < nlocalz + zstart; k++){
		for(int j = ystart; j < nlocaly + ystart; j++){
			for(int i = xstart; i < nlocalx + xstart; i++){
				kk = IndexOf3DIndex( i, j, k);
				if(x[kk * 3] != -100)
					outfile << id[kk] << " "<< x[kk * 3] << " " << x[kk * 3 + 1] << " "<< x[kk * 3 + 2] << std::endl;
			}
		}
	}
	outfile << "print_inter" << std::endl;
	for(int i = 0; i < nlocalinter; i++){
		outfile << idinter[i] << " " << xinter[i][0] << " " << xinter[i][1] << " " << xinter[i][2] << std::endl;
	}
	outfile.close();
}

int atom::getnlocalatom(){
	return (nlocalx * nlocaly * nlocalz);
}

void atom::createphasespace(double factor, int box_x, int box_y, int box_z){
	unsigned long numbefore = 0;
	numbefore += (unsigned long)box_x * box_y * lolocalz;
	numbefore += (unsigned long)lolocaly * box_x * nlocalz;
	numbefore += (unsigned long)lolocalx * nlocaly * nlocalz;
	/*for(int i = 0; i < numbefore; i++){
		uniform();
		uniform();
		uniform();
	}*/
	int xstart = lolocalx - loghostx;
	int ystart = lolocaly - loghosty;
	int zstart = lolocalz - loghostz;
	int kk;
	for(int k = zstart; k < nlocalz + zstart; k++){
		for(int j = ystart; j < nlocaly + ystart; j++){
			for(int i = xstart; i < nlocalx + xstart; i++){
				kk = IndexOf3DIndex( i, j, k);
				id[kk] = ++numbefore;
				kk *= 3;
				x[kk] = (lolocalx + (i - xstart)) * (_latticeconst / 2);
				x[kk + 1] = (lolocaly + (j - ystart)) * _latticeconst + (i % 2) * (_latticeconst / 2);
				x[kk + 2] = (lolocalz + (k - zstart)) * _latticeconst + (i % 2) * (_latticeconst / 2);
				v[kk] = (uniform() - 0.5) * factor;
				v[kk + 1] = (uniform() - 0.5) * factor;
				v[kk + 2] = (uniform() - 0.5) * factor;
			}
		}
	}
}

double atom::uniform(){
	int k = _seed/IQ;
	_seed = IA*(_seed-k*IQ) - IR*k;
	if (_seed < 0) _seed += IM;
	double ans = AM*_seed;
	return ans;
}

void atom::vcm(double mass, double masstotal, double *p){
	double massone;
	int xstart = lolocalx - loghostx;
	int ystart = lolocaly - loghosty;
	int zstart = lolocalz - loghostz;
	int kk;
	for(int k = zstart; k < nlocalz + zstart; k++){
		for(int j = ystart; j < nlocaly + ystart; j++){
			for(int i = xstart; i < nlocalx + xstart; i++){
				kk = IndexOf3DIndex( i, j, k) * 3;
				massone = mass;
				p[0] += v[kk]*massone;
				p[1] += v[kk + 1]*massone;
				p[2] += v[kk + 2]*massone;
			}
		}
	}
}

void atom::zero_momentum(double* vcm){
	int xstart = lolocalx - loghostx;
	int ystart = lolocaly - loghosty;
	int zstart = lolocalz - loghostz;
	int kk;
	for(int k = zstart; k < nlocalz + zstart; k++){
		for(int j = ystart; j < nlocaly + ystart; j++){
			for(int i = xstart; i < nlocalx + xstart; i++){
				kk = IndexOf3DIndex( i, j, k) * 3;
				v[kk] -= vcm[0];
				v[kk + 1] -= vcm[1];
				v[kk + 2] -= vcm[2];
			}
		}
	}
}

double atom::compute_scalar(double mass){
	double t = 0.0;
	int xstart = lolocalx - loghostx;
	int ystart = lolocaly - loghosty;
	int zstart = lolocalz - loghostz;
	int kk;
	for(int k = zstart; k < nlocalz + zstart; k++){
		for(int j = ystart; j < nlocaly + ystart; j++){
			for(int i = xstart; i < nlocalx + xstart; i++){
				kk = IndexOf3DIndex( i, j, k) * 3;
				 t += (v[kk]*v[kk] + v[kk + 1]*v[kk + 1] + v[kk + 2]*v[kk + 2]) * mass;
			}
		}
	}
	return t;
}

void atom::rescale(double scalar, double t_set){
	double factor = sqrt(t_set / scalar);
	int xstart = lolocalx - loghostx;
	int ystart = lolocaly - loghosty;
	int zstart = lolocalz - loghostz;
	int kk;
	for(int k = zstart; k < nlocalz + zstart; k++){
		for(int j = ystart; j < nlocaly + ystart; j++){
			for(int i = xstart; i < nlocalx + xstart; i++){
				kk = IndexOf3DIndex( i, j, k) * 3;
				v[kk] *= factor;
				v[kk + 1] *= factor;
				v[kk + 2] *= factor;
			}
		}
	}
}

