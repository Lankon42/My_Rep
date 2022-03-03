#ifndef H_MY_MATH
#define H_MY_MATH

#include <vector>
#include<iostream>
#include<stdexcept>
#include<string>
#include<cmath>
#include<fstream>
#include <iomanip>
#include <sstream>

class Exception1 : public std::exception{
	public:
		const char* what() const noexcept {
			return "dimensions are inappropriate.";
		}
};

class Exception2 : public std::exception{
	public:
		const char* what() const noexcept {
			return "these matrixes are not vectors.";
		}
};

class Exception3 : public std::exception{
	public:
		const char* what() const noexcept {
			return "one or both vectors have lenght equal to zero.";
		}
};

class Exception4 : public std::exception{
	public:
		const char* what() const noexcept {
			return "square matrix expected.";
		}
};

class Exception5 : public std::exception{
	public:
		const char* what() const noexcept {
			return "invertible matrix expected.";
		}
};

class Exception6 : public std::exception{
	public:
		const char* what() const noexcept {
			return "file does not exist.";
		}
};

class Exception7 : public std::exception{
	public:
		const char* what() const noexcept {
			return "empty file.";
		}
};

class SizeError	// äëÿ ïðîâåðêè ñîîòâåòñòâóþò ëè ðàçìåðû ìàòðèö óñëîâèÿì
{
public:
	std::string size_err;
	SizeError(std::string e) : size_err(e) {};
};

class InputError	// äëÿ ïðîâåðêè ïðàâèëüíîñòè ââîäèìûõ äàííûõ/ïåðåäàâàåìûõ ïàðàìåòðîâ èëè ìàññèâîâ
{
public:
	std::string input_err;
	InputError(std::string e) : input_err(e) {};
};

class ZeroDetError	// äëÿ ïðîâåðêè âûðîæäåííîñòè ìàòðèöû
{
public:
	std::string det_err;
	ZeroDetError(std::string e) : det_err(e) {};
};

class matrix{
	public:
	std::vector<std::vector<double>> matr;
	
	matrix(std::vector<std::vector<double>> vecVec = {{}}); 
	
	virtual matrix mul(double Num);
	virtual matrix Adam(const matrix& M2);
	virtual double det();
	virtual double Fronebius();
	virtual matrix transpose();
	virtual matrix inverse();
	virtual long rank();
	virtual void writeBin(std::ostream& out);
	virtual void readBin(std::istream& in);
	
	double scalar(const matrix& M);
	double norm ();
	double maxnorm ();
	double angle(const matrix& M);
	
	friend std::ostream& operator <<(std::ostream& out, const matrix& M);
	friend std::istream& operator>> (std::istream &in, matrix &M);
	friend matrix operator +(const matrix& M1, const matrix& M2);
	friend matrix operator -(const matrix& M1, const matrix& M2);
	friend matrix operator *(const matrix& M1, const matrix& M2);
};


class identity: public matrix
{
	protected:
	
	public:
	identity(long N = 0);
	
	virtual matrix Adam(const matrix& M2);
	virtual matrix transpose(){return (*this);}
	virtual double det(){return 1;}
	virtual double Fronebius(){return pow(this->matr.size(), 1./2);}
	
	friend std::ostream& operator <<(std::ostream& out, const matrix& Matrix);
	friend std::istream& operator>>(std::istream& in, matrix& M);
	friend matrix operator +(const matrix& M1, const matrix& M2);
	friend matrix operator -(const matrix& M1, const matrix& M2);
	friend matrix operator *(const matrix& M1, const matrix& M2);
};

class diagonal: public matrix
{
	protected:
	
	public:
	diagonal(long N = 0, double* mass = NULL);
	virtual double det();
	virtual matrix transpose(){return (*this);}
	virtual double Fronebius();
	friend std::ostream& operator <<(std::ostream& out, const matrix& Matrix);
	friend std::istream& operator>>(std::istream& in, matrix& M);
	friend matrix operator +(const matrix& M1, const matrix& M2);
	friend matrix operator -(const matrix& M1, const matrix& M2);
	friend matrix operator *(const matrix& M1, const matrix& M2);
};

class upper: public matrix
{
	protected:
	
	public:
	upper(std::vector<std::vector<double>> vecVec = {{}});
	
	virtual double det();
	friend std::ostream& operator <<(std::ostream& out, const matrix& Matrix);
	friend std::istream& operator>>(std::istream& in, matrix& M);
	friend matrix operator +(const matrix& M1, const matrix& M2);
	friend matrix operator -(const matrix& M1, const matrix& M2);
	friend matrix operator *(const matrix& M1, const matrix& M2);
};

class lower: public matrix
{
	protected:
	
	public:
	lower(std::vector<std::vector<double>> vecVec = {{}});
	virtual double det();
	friend std::ostream& operator <<(std::ostream& out, const matrix& Matrix);
	friend std::istream& operator>>(std::istream& in, matrix& M);
	friend matrix operator +(const matrix& M1, const matrix& M2);
	friend matrix operator -(const matrix& M1, const matrix& M2);
	friend matrix operator *(const matrix& M1, const matrix& M2);
};

class symmetric: public matrix
{
	protected:
	
	public:
	symmetric(std::vector<std::vector<double>> vecVec = {{}});
	friend std::ostream& operator <<(std::ostream& out, const matrix& Matrix);
	friend std::istream& operator>>(std::istream& in, matrix& M);
	friend matrix operator +(const matrix& M1, const matrix& M2);
	friend matrix operator -(const matrix& M1, const matrix& M2);
	friend matrix operator *(const matrix& M1, const matrix& M2);
};




#endif // H_MY_MATH

// 1) g++ -c -fPIC my_math.cpp -o my_math.o
// 2) g++ -shared -o libmy_math.so my_math.o
// 3) g++ -c linal.cpp -o linal.o
// 4) g++ linal.o -L/home/emir/Рабочий\ стол/VS/libFolder -lmy_math -o a.out -Wl,-rpath,/home/emir/Рабочий\ стол/VS/libFolder

