#pragma once

#include <random>

class numrand
{
	private:
		int _ir;
		std::mt19937_64 rng_;
		std::uniform_real_distribution<double> dist_;

	public:
		numrand();
		numrand(int ir);
		~numrand();

		double rando();

		void SetIr(int ir);
		int GetIr() const;
};
