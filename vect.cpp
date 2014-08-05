Vect::Vect()
{
}

Vect::Vect(int setN)
{
	N=setN;
	x.resize(N,0);
}

string Vect::typelabel() 
{
	return "vector";
}

Vect Vect::operator%(Vect B)
{
    Vect temp;
    int i;

	if (B.N != N) return (*this);
	
    for (i=0; i<N; i++)
        temp.x[i]=x[i]*B.x[i];

    return temp;
}

Vect Vect::operator+(Vect B)
{
    Vect temp;
    int i;

	if (B.N != N) return (*this);
	
    for (i=0; i<N; i++)
        temp.x[i]=x[i]+B.x[i];

    return temp;
}

Vect Vect::operator-(Vect B)
{
    Vect temp;
    int i;

	if (B.N!=N) return (*this);

    for (i=0; i<N; i++)
        temp.x[i]=x[i]-B.x[i];

    return temp;
}

Vect Vect::operator+=(Vect B)
{
    (*this) = (*this) + B;
    return (*this);
}

Vect Vect::operator-=(Vect B)
{
    (*this) = (*this) - B;
    return (*this);
}

Vect Vect::operator/(Vect B) // Cross product!
{
    Vect temp;
    int i;

	if ((B.N != N)||(N!=3)) return this;

    for (i=0; i<3; i++)
        temp.x[i]=x[(i+1)%3]*B.x[(i+2)%3]-x[(i+2)%3]*B.x[(i+1)%3];

    return temp;
}

float Vect::operator*(Vect B)
{
    float temp=0;
    int i;

	if (B.N != N) return (*this);

    for (i=0; i<N; i++)
        temp+=x[i]*B.x[i];

    return temp;
}

Vect Vect::operator/(float b)
{
    Vect temp;
    int i;

    for (i=0; i<N; i++)
        temp.x[i]=x[i]/b;

    return temp;
}

Vect operator*(float a, Vect B)
{
    Vect temp;
    int i;

    for (i=0; i<N; i++)
        temp.x[i]=B.x[i]*a;

    return temp;
}

Vect Vect::operator*(float b)
{
    Vect temp;
    int i;

    for (i=0; i<3; i++)
        temp.x[i]=x[i]*b;

    return temp;
}

void Vect::Normalize()
{
    float m=sqrt((*this)*(*this));
	if (m<1e-4) m=1e-4;
   
    (*this)=(*this)/m;
}

bool Vect::Within(Vect A, Vect B)
{
    int i;
    bool out=false;

	if ((A.N != N)||(B.N != N)) return false;

    for (i=0; (i<N)&&!out; i++)
    {
        if (x[i]<=A.x[i]) out=true;
        if (x[i]>B.x[i]) out=true;
    }

    return !out;
}

void Vect::zero()
{
	for (int i=0;i<N;i++) x[i]=0;
}

void Vect::resize(int newN)
{
	N=newN;
	x.resize(0,newN);
}
