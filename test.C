#include "rootheader.h"
int hello(int i){

	if(i==3) goto end;
	cout << i <<"\n";
	cout << i <<"\n";

end:
	cout << "the end\n";

}
int main(){

	for(int i=0;i<5;i++){
		hello(i);
	}
	return 0;
}
