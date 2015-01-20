#include "Plane3D.h"
#include "Util.h"
#include <queue>
#include <algorithm>

using namespace std;

namespace Rgbd
{
	bool More_Plane3D(const Plane3D& a,const Plane3D& b)
	{
		return a.inliers > b.inliers;
	}


	Plane3D coff_cloud(pcl::PointCloud<pcl::Normal>::Ptr normals, PointCloud::Ptr cloud)
	{
		//using coff
		double Ex = 0,Ey = 0,Ez = 0,Ex2 = 0,Ey2 = 0,Ez2 = 0,Exy = 0,Eyz = 0,Exz = 0;
		for(int i = 0;i < cloud->size();i++)
		{
			if(!ISNAN(normals->at(i).normal_x))
			{
			
					pcl::PointXYZRGB pp = cloud->at(i);
					//update
					Ex += pp.x;
					Ey += pp.y;
					Ez += pp.z;
					Ex2 += pp.x*pp.x;
					Ey2 += pp.y*pp.y;
					Ez2 += pp.z*pp.z;
					Exy += pp.x*pp.y;
					Eyz += pp.y*pp.z;
					Exz += pp.x*pp.z;

			
			}
		}
		double fem = -Ex2*Ey2*Ez2 + Ex2*Eyz*Eyz + Ez2*Exy*Exy + Ey2*Exz*Exz - 2*Exy*Eyz*Exz;
		Plane3D coff;
		if(fem - 0 < 1e-10 && fem - 0 > -1e-10){
			coff.a = (Exz*Ey2 - Exy*Eyz)/(Exy*Exy - Ex2*Ey2);
			coff.b = -(Exy*Exz - Ex2*Eyz)/(Exy*Exy - Ex2*Ey2);
			coff.c = 1;
			coff.d = 0;
		}else{
			coff.a = (Ex*Eyz*Eyz - Exz*Ey*Eyz - Ex*Ey2*Ez2 + Exy*Ey*Ez2 + Exz*Ey2*Ez - Exy*Eyz*Ez)/fem;
			coff.b = (Exz*Exz*Ey - Ex*Exz*Eyz + Ex*Exy*Ez2 - Exy*Exz*Ez - Ex2*Ey*Ez2 + Ex2*Eyz*Ez)/fem;
			coff.c = (Exy*Exy*Ez + Ex*Exz*Ey2 - Ex*Exy*Eyz - Exy*Exz*Ey + Ex2*Ey*Eyz - Ex2*Ey2*Ez)/fem;
			coff.d = 1;
		}
		double u = NORM2(coff.a, coff.b, coff.c);
		coff.a /= u;
		coff.b /= u;
		coff.c /= u;
		coff.d /= u;

		//cout<<"("<<coff.a<<","<<coff.b<<","<<coff.c<<","<<coff.d<<") "<<"inliers:"<<coff.inliers<<endl;
		return coff;
	}

	void dbplane_segment(PointCloud::Ptr cloud, PointNormal::Ptr normals,
		                 int minInlier, double sthre, double athre, double dthre,
						 std::vector<Plane3D>& coffs, vector<int>& seg)
	{
		//int minInlier = 1000;
		vector<int> valid;		
		int csize = 0;

		for(int i = 0;i < cloud->size();i++)
		{
			if(seg[i] == -1 && !ISNAN(normals->at(i).normal_x))
			{
				int id = csize;
				int count = 1;
				seg[i] = id;
				queue<int> search;
				search.push(i);

				double Ex = 0,Ey = 0,Ez = 0,Ex2 = 0,Ey2 = 0,Ez2 = 0,Exy = 0,Eyz = 0,Exz = 0;

				double avern[3];
				int averc = 0;
				avern[0] = 0;
				avern[1] = 0;
				avern[2] = 0;
				while(!search.empty())
				{
					int p = search.front();
					search.pop();
					int y = p / cloud->width;
					int x = p - y*cloud->width;
					pcl::Normal np = normals->at(p);
					avern[0] = (avern[0]*averc + np.normal_x)/(averc + 1);
					avern[1] = (avern[1]*averc + np.normal_y)/(averc + 1);
					avern[2] = (avern[2]*averc + np.normal_z)/(averc + 1);
					averc++;
					double norm = NORM2(avern[0], avern[1], avern[2]);
					avern[0] /= norm;
					avern[1] /= norm;
					avern[2] /= norm;

					pcl::PointXYZRGB pp = cloud->at(p);


					Ex += pp.x;
					Ey += pp.y;
					Ez += pp.z;
					Ex2 += pp.x*pp.x;
					Ey2 += pp.y*pp.y;
					Ez2 += pp.z*pp.z;
					Exy += pp.x*pp.y;
					Eyz += pp.y*pp.z;
					Exz += pp.x*pp.z;

					for(int jx = x-1;jx <= x+1;jx++)
					{
						if(jx < 0) continue;
						if(jx >= cloud->width) continue;
						for(int jy = y-1;jy <= y+1;jy++)
						{
							if(jy < 0) continue;
							if(jy >= cloud->height) continue;
							int nid = jy*cloud->width+jx;
							if(seg[nid] == -1 && !ISNAN(normals->at(nid).normal_x))
							{
								pcl::Normal nn = normals->at(nid);
								pcl::PointXYZRGB pn = cloud->at(nid);
								//double disn = sqrt((np.normal_x-nn.normal_x)*(np.normal_x-nn.normal_x)+(np.normal_y-nn.normal_y)*(np.normal_y-nn.normal_y)+(np.normal_z-nn.normal_z)*(np.normal_z-nn.normal_z));
								double angn = avern[0] * nn.normal_x + avern[1] * nn.normal_y + avern[2] * nn.normal_z;
								double disp = NORM2(pp.x-pn.x, pp.y-pn.y, pp.z-pn.z);
								//double dis = sqrt(disn+disp);
								if(angn > athre && disp < sthre)
								{
									search.push(nid);
									seg[nid] = id;
									count++;
								}
							}
						}
					}
				}
				double fem = -Ex2*Ey2*Ez2 + Ex2*Eyz*Eyz + Ez2*Exy*Exy + Ey2*Exz*Exz - 2*Exy*Eyz*Exz;
				Plane3D coff;
				if(fem - 0 < 1e-10 && fem - 0 > -1e-10){
					coff.a = (Exz*Ey2 - Exy*Eyz)/(Exy*Exy - Ex2*Ey2);
					coff.b = -(Exy*Exz - Ex2*Eyz)/(Exy*Exy - Ex2*Ey2);
					coff.c = 1;
					coff.d = 0;
				}else{
					coff.a = (Ex*Eyz*Eyz - Exz*Ey*Eyz - Ex*Ey2*Ez2 + Exy*Ey*Ez2 + Exz*Ey2*Ez - Exy*Eyz*Ez)/fem;
					coff.b = (Exz*Exz*Ey - Ex*Exz*Eyz + Ex*Exy*Ez2 - Exy*Exz*Ez - Ex2*Ey*Ez2 + Ex2*Eyz*Ez)/fem;
					coff.c = (Exy*Exy*Ez + Ex*Exz*Ey2 - Ex*Exy*Eyz - Exy*Exz*Ey + Ex2*Ey*Eyz - Ex2*Ey2*Ez)/fem;
					coff.d = 1;
				}
				double u = NORM2(coff.a, coff.b, coff.c);
				coff.a /= u;
				coff.b /= u;
				coff.c /= u;
				coff.d /= u;
				/*
				if(count > minInlier)
				{
					for(int j = 0;j < cloud->size();j++)
					{
						if((seg[j] == -1 || (seg[j] < valid.size() && valid[seg[j]] == -1)) && !isnan(normals->at(j).normal_x))
						{
							pcl::PointXYZRGB pj = cloud->at(j);
							if(abs(coff.a * pj.x + coff.b * pj.y + coff.c * pj.z - coff.d) < dthre)
							{
								seg[j] = id;
							
								Ex += pj.x;
								Ey += pj.y;
								Ez += pj.z;
								Ex2 += pj.x*pj.x;
								Ey2 += pj.y*pj.y;
								Ez2 += pj.z*pj.z;
								Exy += pj.x*pj.y;
								Eyz += pj.y*pj.z;
								Exz += pj.x*pj.z;
								count++;
							}
						}
					}
					double fem2 = -Ex2*Ey2*Ez2 + Ex2*Eyz*Eyz + Ez2*Exy*Exy + Ey2*Exz*Exz - 2*Exy*Eyz*Exz;
					if(fem2 - 0 < 1e-10 && fem2 - 0 > -1e-10){
						coff.a = (Exz*Ey2 - Exy*Eyz)/(Exy*Exy - Ex2*Ey2);
						coff.b = -(Exy*Exz - Ex2*Eyz)/(Exy*Exy - Ex2*Ey2);
						coff.c = 1;
						coff.d = 0;
					}else{
						coff.a = (Ex*Eyz*Eyz - Exz*Ey*Eyz - Ex*Ey2*Ez2 + Exy*Ey*Ez2 + Exz*Ey2*Ez - Exy*Eyz*Ez)/fem;
						coff.b = (Exz*Exz*Ey - Ex*Exz*Eyz + Ex*Exy*Ez2 - Exy*Exz*Ez - Ex2*Ey*Ez2 + Ex2*Eyz*Ez)/fem;
						coff.c = (Exy*Exy*Ez + Ex*Exz*Ey2 - Ex*Exy*Eyz - Exy*Exz*Ey + Ex2*Ey*Eyz - Ex2*Ey2*Ez)/fem;
						coff.d = 1;
					}
					u = sqrt(coff.a*coff.a + coff.b*coff.b + coff.c*coff.c);
					coff.a /= u;
					coff.b /= u;
					coff.c /= u;
					coff.d /= u;
				}
				*/
				coff.inliers = count;
				coff.id = id;
				csize++;
				if(coff.inliers > minInlier)
				{
					
					coffs.push_back(coff);
					valid.push_back(coffs.size()-1);
				}
				else
				{
					valid.push_back(-1);
				}
			}
		}
		sort(coffs.begin(), coffs.end(),More_Plane3D);
		for(int i = 0;i < coffs.size();i++)
		{
			valid[coffs[i].id] = i;
		}
		for(int i = 0;i < cloud->size();i++)
		{
			if(seg[i] >= 0) seg[i] = valid[seg[i]];
			
		}
		
		for(int i =0;i < coffs.size();i++) cout<<"("<<coffs[i].a<<","<<coffs[i].b<<","<<coffs[i].c<<","<<coffs[i].d<<") "<<"inliers:"<<coffs[i].inliers<<endl;
		cout<<endl;
	}

	Plane3D coffCal::getCoff()	
	{
		double fem = -Ex2*Ey2*Ez2 + Ex2*Eyz*Eyz + Ez2*Exy*Exy + Ey2*Exz*Exz - 2*Exy*Eyz*Exz;
		Plane3D coff;
		if(fem - 0 < 1e-10 && fem - 0 > -1e-10){
			coff.a = (Exz*Ey2 - Exy*Eyz)/(Exy*Exy - Ex2*Ey2);
			coff.b = -(Exy*Exz - Ex2*Eyz)/(Exy*Exy - Ex2*Ey2);
			coff.c = 1;
			coff.d = 0;
		}else{
			coff.a = (Ex*Eyz*Eyz - Exz*Ey*Eyz - Ex*Ey2*Ez2 + Exy*Ey*Ez2 + Exz*Ey2*Ez - Exy*Eyz*Ez)/fem;
			coff.b = (Exz*Exz*Ey - Ex*Exz*Eyz + Ex*Exy*Ez2 - Exy*Exz*Ez - Ex2*Ey*Ez2 + Ex2*Eyz*Ez)/fem;
			coff.c = (Exy*Exy*Ez + Ex*Exz*Ey2 - Ex*Exy*Eyz - Exy*Exz*Ey + Ex2*Ey*Eyz - Ex2*Ey2*Ez)/fem;
			coff.d = 1;
		}
		double u = NORM2(coff.a, coff.b, coff.c);
		coff.a /= u;
		coff.b /= u;
		coff.c /= u;
		coff.d /= u;

		return coff;
	}

	double p_distance(Plane3D& a, Plane3D& b)
	{
		return sqrt((a.a - b.a) * (a.a - b.a) + (a.b - b.b) * (a.b - b.b) + (a.c - b.c) * (a.c - b.c) + (a.d - b.d) * (a.d - b.d));
	}

	double p_distance(Plane3D& p, double a, double b, double c, double d)
	{
		return sqrt((p.a - a) * (p.a - a) + (p.b - b) * (p.b - b) + (p.c - c) * (p.c - c) + (p.d - d) * (p.d - d));
	}

	double p_angle(Plane3D& a, Plane3D& b)
	{
		return (a.a * b.a + a.b * b.b + a.c * b.c)/sqrt(a.a*a.a + a.b*a.b + a.c*a.c)/sqrt(b.a*b.a + b.b*b.b + b.c*b.c);
	}

	double p_angle(Plane3D& p, double a, double b, double c, double d)
	{
		return (p.a * a + p.b * b + p.c * c)/sqrt(p.a*p.a + p.b*p.b + p.c*p.c)/sqrt(a*a + b*b + c*c);
	}
}