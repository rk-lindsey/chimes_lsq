#include<iostream>
#include<vector>

using namespace std;

int main()
{

	vector<int>myvec(8);
	cout << myvec.size() << endl;
	vector<int>().swap(myvec);
	cout << myvec.size() << endl;
}
