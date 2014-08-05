class Vect: public DataType
{
	public:
		int N;
		vector<float> x;
		
		string typelabel();		
		Vect();		
		Vect(int setN);
		Vect operator%(Vect B);
		Vect operator+(Vect B);
		Vect operator-(Vect B);
		Vect operator+=(Vect B);
		Vect operator-=(Vect B);
		Vect operator/(Vect B); // Cross product!
		float operator*(Vect B);
		Vect operator/(float b);
		Vect operator*(float b);
		void Normalize();
		bool Within(Vect A, Vect B);
		void zero();
		void resize(int newN);
};

extern Vect operator*(float a, Vect B);
