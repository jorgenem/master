#include <iostream>
using namespace std;
#include <vector>
int main() {
	vector <int> balle = {0,1,2,3,4,5,6,7,8,9,10};
	for (auto it = balle.begin()+2; it != balle.end(); it++) {
		cout << *it << endl;
	}
	return 0;
}