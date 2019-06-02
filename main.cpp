#include "fixed.hpp"
#include "bmp.hpp"
const fixed<1,7> center = {0,0,0,0,0,0,0,0};
const fixed<1,7> vertical = {1,0,0,0,0,0,0,0};
const fixed<1,7> horizontal = {1,0,0,0,0,0,0,0};
int main(){
	fixed<1, 1> fl(0.43d);
	std::cout << fl.to_gmp_float() << "\n";
	/*fixed<1, 7> a;
	std::vector<vec3> colors(1000 * 1000);
	for(size_t i = 0;i < 1000;i++){
		for(size_t j = 0;j < 1000;j++){
			
		}
	}*/
}
