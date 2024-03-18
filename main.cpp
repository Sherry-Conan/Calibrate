#include"Eu.h"

Int_t main(){

	std::ofstream GeFile("kbGe3092.txt"); // 打开文件
	std::ofstream LaBrFile("LaBr3092.txt"); // 打开文件

    Eu *ana = new Eu(3092); 

	ana -> GetBAndKmain(); 

	ana -> WriteOut(GeFile, LaBrFile); 

	return 0; 
}