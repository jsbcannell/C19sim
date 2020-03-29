/* c19sim.cpp
 * Copyright (C) 2017-2019  Jake S. Cannell
*/

/* GNU Affero General Public License, Version 3 {{{ */
/*
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.

 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.

 * You should have received a copy of the GNU Affero General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/




#include <string>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <sstream>
#include <random>
#include <stdarg.h>
#include <time.h>
#include <chrono>
#include <random>

using namespace std;

class CSVRow
{
    public:
        std::string const& operator[](std::size_t index) const
        {
            return m_data[index];
        }
        std::size_t size() const
        {
            return m_data.size();
        }
        void readNextRow(std::istream& str)
        {
            std::string         line;
            std::getline(str, line);

            std::stringstream   lineStream(line);
            std::string         cell;

            m_data.clear();
            while(std::getline(lineStream, cell, ','))
            {
                m_data.push_back(cell);
            }
            // This checks for a trailing comma with no data after it.
            if (!lineStream && cell.empty())
            {
                // If there was a trailing comma then add an empty element.
                m_data.push_back("");
            }
        }
    private:
        std::vector<std::string>    m_data;
};


std::istream& operator>>(std::istream& str, CSVRow& data)
{
    data.readNextRow(str);
    return str;
}

template <class GT, class DeathTimeDistT>
int sim_run(GT& G, const vector<int>& cases, double IFR, DeathTimeDistT dtd, int today, int obs_deaths)
{
		int sim_deaths(0);
		uniform_real_distribution<> uni_dist(0.0, 1.0);
		for (int i(0); i < cases.size(); i++)
		{
			double total_infected = cases[i];
			//printf("total_infected[%i]: %f \n", i, total_infected);

			for (int j(0); j < total_infected; j++)
			{
				//int infection_time = i - round( ctd(G) );
				if (uni_dist(G) < IFR) { // death check
					int death_time = i + round( dtd(G) );
					//printf("%i ", death_time);
					if (death_time <= today) {
						++sim_deaths;
						//if (++sim_deaths > obs_deaths) return false;
					}
				}
			} // j
			//printf("\n");
		} // i
		return sim_deaths;
} // sim_run

template <class GT, class DeathTimeDistT>
double run_sims(GT& G, const vector<int>& cases, double IFR, DeathTimeDistT dtd, int today, int obs_deaths, int N)
{
	int nmatches(0);
	for (int i(0); i < N; i++) {
		int sim_deaths = sim_run(G,cases,IFR,dtd,today,obs_deaths);
		if (sim_deaths == obs_deaths) nmatches++;
	}
	return double(nmatches) / double(N);
}


void simulate(double IFR, double C_N, int obs_deaths, int shape)
{
	std::ifstream       file("iceland_cases.csv");

	vector<int> cases;

	int cases_NUHI(0);
	int cases_DGEN(0);

	int i_(0);
	CSVRow              row;
	while(file >> row)
	{
			//printf("%s %s %s \n", row[0].c_str(), row[1].c_str(), row[2].c_str());
			if (i_ > 0) {
				int n0(0), n1(0);
				if (row[1].size() > 2) n0 = stoi( row[1].substr(1,row[1].size()-1) );
				if (row[2].size() > 2) n1 = stoi( row[2].substr(1,row[2].size()-1) );
				//printf("%i %i \n", n0, n1);
				cases.push_back(n0 + n1);
				cases_NUHI += n0;
				cases_DGEN += n1;
			}
			i_++;
	}
	int cases_tot = cases_NUHI + cases_DGEN;
	printf("today: %i total cases %i = %i + %i \n", int(cases.size()), cases_tot, cases_NUHI, cases_DGEN);

	double extra_infected = cases_tot * (C_N - 1.0);

	vector<double> extra_cases;

	double ytot = 0.0;
	for (int i(0); i < int(cases.size()); i++)
	{
		double x = double(i);
		double sx = 0.3*x -1.0;
		double y = exp(sx) / (exp(sx) + 1.0);
		//printf(" %f", y);
		extra_cases.push_back(y);
		ytot += y;
	}
	//printf("\n");

	for (int i(0); i < int(extra_cases.size()); i++)
	{
		extra_cases[i] = round( extra_cases[i] / ytot * extra_infected );
	}

	for (int i(0); i < cases.size(); i++)
	{
		if (shape == 0) cases[i] += extra_infected / double(cases.size());
		if (shape == 1) cases[i] += int( extra_cases[i] );
		if (shape == 2) cases[i] = int(double(cases[i]) * C_N);
	}

	random_device rd;  //Will be used to obtain a seed for the random number engine
	mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()


	double pval = run_sims(gen, cases, IFR, lognormal_distribution<>(2.8,0.45), cases.size(), obs_deaths, 10000);

	printf("%f \n", pval);

}


// ============================ ===================================

int main(int argc, char* argv[])
{
	printf("%s\n", argv[1]);
	int ND  	 = stoi(argv[1]);
	double IFR = stof(argv[2]);
	double C_N = stof(argv[3]);
	int shape  = stoi(argv[4]);

	printf("c19sim IFR=%f C_N=%f ND=%i shape=%i \n", IFR, C_N, ND, shape);

	simulate(IFR, C_N, ND, shape);


	printf("\n");

	printf("Done!\n");

}


/*

Age  		C  	 N  	C/N
 0- 9		21  50K  0.4/K
10-19		76  45K  1.6/K
20-29  158  50K  3.1/K
30-39	 171  48K  3.5/K
40-49  226  45K  5.0/K
50-59  190  46K  4.1/K
60-69  131  35K  3.7/K
70-79   37  22K  1.6/K
80-89    4  13K  0.3/K
90+      6   3K  2.0/K




*/
