#pragma once

class numrand
{
	private:
		int _ir;

	public:
		numrand();
		numrand(int ir);
		~numrand();

		double rando();

		void SetIr(int ir);
		int GetIr() const;
};
