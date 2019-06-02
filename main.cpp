#include "fixed.hpp"
#include <iostream>
#include <vector>
#include <numeric>
#include <iomanip>
void superfunction(fixed<20,20>& a, fixed<20,20>& b){
	a += b;
}
uint64_t nanoTime(){
	using namespace std;
	using namespace std::chrono;
	return duration_cast<nanoseconds>(high_resolution_clock::now().time_since_epoch()).count();
}
int main(){
	xoshiro_256 gen(42);
	std::cout << std::setprecision(1000);
	fixed<1, 2> a;
	std::generate(a.bits.begin(), a.bits.end(), std::ref(gen));
	
	std::vector<uint64_t> times;
	times.reserve(2002);
	uint64_t mintime = -1;
	for(;;){
		auto t1 = nanoTime();
		for(int i = 0;i < 1000;i++){
			a = a * a;
		}
		auto t2 = nanoTime();
		if((t2 - t1) < mintime){
			mintime = (t2 - t1);
			times.clear();
			times.push_back(t2 - t1);
		}
		else{
			times.push_back(t2 - t1);
			if(times.size() >= 2000){
				break;
			}
		}
	}
	uint64_t totaltime = std::accumulate(times.begin(), times.end(), uint64_t(0));
	std::cout << "1000 Multiplications:\n";
	std::cout << "Min: " << mintime << " ns, Average: " << totaltime / 2000.0 << " ns\n";
	//std::cout << (t2 - t1) / 1000 / 1000.0 << " ms\n" << (a.bits[0] & 1) << "\n";
	return 0;
}
