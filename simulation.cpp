#include <iostream>
using namespace std;

template <class T = int>
class simulation
{
private:
	T dx, dy;
	T x,y,h;
public:
	simulation(T DX=10,T DY=10,T X=500,T Y=500,T H=10)
	{
		dx = DX; dy = DY; x = X; y = Y; h = H;
		cout << x << endl;
	}
	
	void setDepth(T NEWH)
	{
		h = NEWH;
	}

	//show(){}
	
	void set()
	{
		int x;
		T temp;
		cout << "What do you want to set?" << endl;
		cout << "1.dx	2.dy	3.x	4.y	5.h" << endl;
		cin >> x;
		if (x==1) 
		{
			cout << "dx=? " << endl;
			cin >> temp;
			dx = temp;	
		}

		else if (x==2) 
                {
                        cout << "dy=? " << endl;
                        cin >> temp;
                        dy = temp;
                }

		else if (x==3) 
                {
                        cout << "x=? " << endl;
                        cin >> temp;
                        x = temp;
                }

		else if (x==4) 
                {
                        cout << "y=? " << endl;
                        cin >> temp;
                        y = temp;
                }

		else if (x==5) 
                {
                        cout << "h=? " << endl;
                        cin >> temp;
                        h = temp;
                }
		
		else
		{
			cout << "Your entrance is out of list" << endl;
		}	

		
	}

	void get()
	{
		cout << "-----Parameters-----" << endl;
		cout << "dx= " << dx << endl;
		cout << "dy= " << dy << endl;
		cout << "x= " << x  << endl;
		cout << "y= " << y << endl;
		cout << "h= " << h << endl;
	}
	
};

int main()
{
	simulation<>v1;
	v1.get();
	v1.set();
	v1.get();

	return 0;
}
