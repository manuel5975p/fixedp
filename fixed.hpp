#ifndef FIXED_HPP
#define FIXED_HPP
#include <cstddef>
#include <cstdint>
#include <array>
#include <climits>
#include <x86intrin.h>
#include <immintrin.h>
#include <algorithm>
#include <deque>
#include <string>
#include <sstream>
#ifdef GMP_CONVERSIONS
#include <gmpxx.h>
#endif
#include <iostream>
#include <functional>
using std::size_t;
using std::uint64_t;
using uint128_t = unsigned __int128;
using int128_t = __int128;
template<typename T>
struct larger_int_struct{};
template<> struct larger_int_struct<std::uint8_t >{using type = std::uint16_t; };
template<> struct larger_int_struct<std::uint16_t>{using type = std::uint32_t; };
template<> struct larger_int_struct<std::uint32_t>{using type = std::uint64_t; };
template<> struct larger_int_struct<std::uint64_t>{using type =      uint128_t;};
template<> struct larger_int_struct<std::int8_t  >{using type = std::int16_t;  };
template<> struct larger_int_struct<std::int16_t >{using type = std::int32_t;  };
template<> struct larger_int_struct<std::int32_t >{using type = std::int64_t;  };
template<> struct larger_int_struct<std::int64_t >{using type =      int128_t; };
template<typename T>
using larger_int = typename larger_int_struct<T>::type;
template<typename int_t> char adc_intrinsic(char in_flag, int_t a, int_t b, int_t* stor){return 0;};

template<>
char adc_intrinsic<uint64_t>(char in_flag, uint64_t a, uint64_t b, uint64_t* stor){
	return _addcarry_u64(in_flag, a, b, (unsigned long long*)stor);
};
template<>
char adc_intrinsic<uint32_t>(char in_flag, uint32_t a, uint32_t b, uint32_t* stor){
	return _addcarry_u32(in_flag, a, b, (unsigned int*)stor);
};
template<typename int_t> char sbb_intrinsic(char in_flag, int_t a, int_t b, int_t* stor){return 0;};

template<>
char sbb_intrinsic<uint64_t>(char in_flag, uint64_t a, uint64_t b, uint64_t* stor){
	return _subborrow_u64(in_flag, a, b, (unsigned long long*)stor);
};
template<>
char sbb_intrinsic<uint32_t>(char in_flag, uint32_t a, uint32_t b, uint32_t* stor){
	return _subborrow_u32(in_flag, a, b, (unsigned int*)stor);
};

template<typename int_t> int_t mul_intrinsic(int_t a, int_t b, int_t* hi){return 0;};

template<>
uint64_t mul_intrinsic<uint64_t>(uint64_t a, uint64_t b, uint64_t* hi){
	return _mulx_u64(a, b, (unsigned long long*)hi);
};
template<>
uint32_t mul_intrinsic<uint32_t>(uint32_t a, uint32_t b, uint32_t* hi){
	uint64_t _a = a;
	_a *= (uint64_t)b;
	*hi = _a >> 32;
	return _a;
};
template<typename T>
struct float_info{
};
template<>
struct float_info<float>{
	std::uint32_t mantissa;
	std::int16_t exponent;
	bool sign;
	float_info(float a){
		std::uint32_t val = *reinterpret_cast<std::uint32_t*>(&a);
		mantissa = val & ((1u << 23) - 1);
		val >>= 23;
		exponent = val & ((1u << 8) - 1);
		val >>= 8;
		sign = !!val;
	}
};
template<>
struct float_info<double>{
	std::uint64_t mantissa;
	std::int16_t exponent;
	bool sign;
	float_info(double a){
		std::uint64_t val = *reinterpret_cast<std::uint64_t*>(&a);
		mantissa = val & ((1ULL << 52) - 1);
		val >>= 52;
		exponent = val & ((1ULL << 11) - 1);
		val >>= 11;
		sign = !!val;
	}
};
template<typename int_t>
int_t lzcnt_intrinsic(int_t x){

}
template<typename int_t>
int_t tzcnt_intrinsic(int_t x){

}
template<typename int_t>
int_t popcnt_intrinsic(int_t x){

}
template<>
uint64_t lzcnt_intrinsic<uint64_t>(uint64_t x){
	return _lzcnt_u64(x);
}
template<>
uint64_t tzcnt_intrinsic<uint64_t>(uint64_t x){
	return _tzcnt_u64(x);
}
template<>
uint64_t popcnt_intrinsic<uint64_t>(uint64_t x){
	return _popcnt64(x);
}
template<>
uint32_t lzcnt_intrinsic<uint32_t>(uint32_t x){
	return _lzcnt_u32(x);
}
template<>
uint32_t tzcnt_intrinsic<uint32_t>(uint32_t x){
	return _tzcnt_u32(x);
}
template<>
uint32_t popcnt_intrinsic<uint32_t>(uint32_t x){
	return _popcnt32(x);
}
template<size_t n, typename int_t>
struct pointer_adder{
	pointer_adder(){}
	void operator()(int_t* first, const int_t* second, bool carry){
		pointer_adder<n - 1, int_t>()(first, second, adc_intrinsic<int_t>(carry, first[n], second[n], first + n));
	}
};

template<typename int_t>
struct pointer_adder<0, int_t>{
	pointer_adder(){}
	void operator()(int_t* first, const int_t* second, bool carry){
		adc_intrinsic<int_t>(carry, first[0], second[0], first);
	}
};

template<size_t n, typename int_t>
struct pointer_subtractor{
	pointer_subtractor(){}
	void operator()(int_t* first, const int_t* second, bool carry){
		pointer_subtractor<n - 1, int_t>()(first, second, sbb_intrinsic<int_t>(carry, first[n], second[n], first + n));
	}
};

template<typename int_t>
struct pointer_subtractor<0, int_t>{
	pointer_subtractor(){}
	void operator()(int_t* first, const int_t* second, bool carry){
		sbb_intrinsic<int_t>(carry, first[0], second[0], first);
	}
};

template<size_t m, size_t n, typename int_t = uint64_t>
struct fixed{
	std::array<int_t, m + n> bits;
	static constexpr size_t M = m;
	static constexpr size_t N = n;
	fixed(std::initializer_list<int_t> init) : fixed(){
		std::copy(init.begin(), init.end(), bits.begin());
	}
	fixed(){
		std::fill(bits.begin(), bits.end(), 0);
	}
	fixed(double x) : fixed(){
		float_info<double> info(x);
		std::uint64_t val = *reinterpret_cast<std::uint64_t*>(&x);
		
		bits[m] = info.mantissa << 12;
		bits[m - 1] = 1;
		shiftLeft(info.exponent);
	}
	template<typename RNG>
	fixed(RNG& gen, bool dummy){
		std::generate(bits.begin(), bits.end(), std::ref(gen));
	}
	fixed& operator+=(const fixed<m,n,int_t>& o){
		pointer_adder<m + n - 1, int_t>()(bits.data(), o.bits.data(), 0);
		return *this;
	}
	fixed& operator-=(const fixed<m,n,int_t>& o){
		pointer_subtractor<m + n - 1, int_t>()(bits.data(), o.bits.data(), 0);
		return *this;
	}
	
	fixed& csleft(int _n){
		if(_n < 0)return csright(-_n);
		if(n == 0)return *this;
		if(_n >= m + n){
			std::fill(bits.begin(), bits.end(), 0);
			return *this;
		}
		for(size_t i = 0;i < m + n - _n;i++){
			bits[i] = bits[i + _n];
		}
		for(size_t i = m + n - _n;i < m + n;i++){
			bits[i] = 0;
		}
		return *this;
	}
	fixed& csright(int _n){
		if(_n < 0)return csleft(-_n);
		if(n == 0)return *this;
		if(_n >= m + n){
			std::fill(bits.begin(), bits.end(), 0);
			return *this;
		}
		for(size_t i = m + n - 1;i >= _n;i--){
			bits[i] = bits[i - _n];
		}
		for(size_t i = 0;i < _n;i++){
			bits[i] = 0;
		}
		return *this;
	}
	fixed& operator*=(int_t mul){
		mult_comp<m + n - 1>(mul, 0);
		return *this;
	}
	
	fixed operator*(const fixed<m, n, int_t>& o)const{
		fixed<m, n> accum;
		for(ssize_t offset = 0;offset < m;offset++){
			int_t carry = 0;
			for(ssize_t i = m + n - 1;i >= offset;i--){
				int_t hi = 0;
				int_t prod = mul_intrinsic(o.bits[m - offset - 1], bits[i], &hi);
				char add_overflow = !!adc_intrinsic(0, prod, carry, &prod);
				add_overflow += !!adc_intrinsic(0, prod, accum.bits[i - offset], &accum.bits[i - offset]);
				carry = hi + add_overflow;
			}
		}
		
		for(ssize_t offset = 0;offset < n;offset++){
			int_t carry;
			mul_intrinsic(o.bits[m + offset], bits[m + n - offset - 1], &carry);
			for(ssize_t i = m + n - offset - 2;i >= 0;i--){
				int_t hi;
				int_t prod = mul_intrinsic(o.bits[m + offset], bits[i], &hi);
				char add_overflow = !!adc_intrinsic(0, prod, carry, &prod);
				add_overflow += !!adc_intrinsic(0, prod, accum.bits[i + offset + 1], &accum.bits[i + offset + 1]);
				carry = hi + add_overflow;
			}
			size_t o = offset;
			while((carry = adc_intrinsic(0, accum.bits[o], carry, &accum.bits[o]))){--o;}
		}
		return accum;
	}

	fixed operator/(const fixed<m, n>& o)const{
		int64_t shift = 0;
		fixed<m, n> subber(o);
		fixed<m, n> tcopy(*this);
		fixed<m, n> result;
		
		while(!tcopy.isZero()){
			int64_t cshift = int64_t(int64_t(subber.lzcnt()) - int64_t(tcopy.lzcnt()));
			subber <<= cshift;
			if(subber > tcopy){
				subber >>= 1;
				cshift--;
			}
			shift += cshift;
			if(-shift > int64_t(n) * int64_t(sizeof(int_t)) * int64_t(CHAR_BIT)){
				break;
			}
			tcopy -= subber;
			result.putBit(shift);
		}
		return result;
	}
	bool isZero()const{
		for(auto it = rbegin();it != rend();it++){
			if(*it)return false;
		}
		return true;
	}
	template<size_t i>
	void mult_comp(int_t mul, int_t carry){
		int_t _carry = 0;
		bits[i] = mul_intrinsic(bits[i], mul, &_carry);
		bits[i] += carry;
		if(bits[i] < carry)
			_carry++;
		if constexpr(i > 0){
			mult_comp<i - 1>(mul, _carry);
		}
	}
	void putBit(int64_t offset){
		int64_t chunk = (offset - int64_t(sizeof(int_t)) * int64_t(CHAR_BIT) + 1) / int64_t(sizeof(int_t) * CHAR_BIT);
		bits[int64_t(m) - 1 - (offset - int64_t(sizeof(int_t)) * int64_t(CHAR_BIT) + 1) / int64_t(sizeof(int_t) * CHAR_BIT)] |= (int_t(1) << (offset % (sizeof(int_t) * CHAR_BIT)));
	}
	int_t operator%(int_t mod)const{
		return smallmod<0>(mod, 0);
	}
	
	template<size_t i>
	int_t smallmod(int_t mod, larger_int<int_t> carry)const{
		if constexpr(i < m - 1){
			return smallmod<i + 1>(mod, (larger_int<int_t>)(((larger_int<int_t>)bits[i] + carry) % mod) << sizeof(int_t) * CHAR_BIT);
		}
		return ((larger_int<int_t>)bits[m - 1] + carry) % mod;
	}
	
	fixed& operator/=(int_t div){
		smalldiv<0>(div, 0);
		return *this;
	}
	
	template<size_t i>
	void smalldiv(int_t div, larger_int<int_t> carry){
		carry += bits[i];
		larger_int<int_t> ncarry = carry % div;
		bits[i] = carry / div;
		if constexpr(i < m - 1){
			smalldiv<i + 1>(div, ncarry << sizeof(int_t) * CHAR_BIT);
		}
	}
	fixed operator/=(const fixed<m,n, int_t>& o)const{
		auto x = *this / o;
		std::copy(x.begin(), x.end(), begin());
		return this;
	}
	fixed operator/(int_t o)const{
		fixed<m,n,int_t> copy(*this);
		copy /= o;
		return copy;
	}
	fixed operator+(const fixed<m,n, int_t>& o)const{
		fixed<m,n,int_t> copy(*this);
		copy += o;
		return copy;
	}
	fixed operator-(const fixed<m,n, int_t>& o)const{
		fixed<m,n,int_t> copy(*this);
		copy -= o;
		return copy;
	}
	fixed& operator*=(const fixed<m,n, int_t>& o){
		auto x = *this * o;
		std::copy(x.begin(), x.end(), begin());
		return this;
	}
	#ifdef GMP_CONVERSIONS
	mpf_class to_gmp_float()const{
		mpf_class ret(0, bits.size() * sizeof(int_t) * CHAR_BIT + 64);
		for(ssize_t i = 0;i < bits.size();i++){
			mpf_class two(2, bits.size() * sizeof(int_t) * CHAR_BIT + 64);
			mpf_class one(1, bits.size() * sizeof(int_t) * CHAR_BIT + 64); 
			mpf_pow_ui(two.get_mpf_t(), two.get_mpf_t(), sizeof(int_t) * CHAR_BIT * std::abs(ssize_t(m) - i - 1));
			if(ssize_t(m) - i - 1 < 0){
				mpf_div(two.get_mpf_t(), one.get_mpf_t(), two.get_mpf_t());
			}
			two *= bits[i];
			ret += two;
		}
		return ret;
	}
	#endif
	template<size_t i = m - 1>
	bool has_bits_before()const{
		if(bits[i]){
			return true;
		}
		if constexpr(i > 0){
			return has_bits_before<i - 1>();
		}
		return false;
	}
	
	template<size_t i = m>
	bool has_bits_after()const{
		if(bits[i])return true;
		if constexpr(i < m + n - 1){
			return has_bits_before<i + 1>();
		}
		return false;
	}
	std::string to_expression_string()const{
		int exp = (m - 1) * 64;
		std::stringstream str;
		for(size_t i = 0;i < m + n;i++){
			str << bits[i] << " * 2^" << exp;
			exp -= 64;
			if(i < m + n - 1)
				str << " + ";
		}
		return str.str();
	}
	std::string to_string()const{
		std::deque<char> chars;
		auto num = *this;
		while(num.has_bits_before()){
			chars.push_front('0' + (num % 10));
			num /= 10;
		}
		if(chars.empty()){
			chars.push_back('0');
		}
		if(has_bits_after()){
			chars.push_back('.');
			num = *this;
			for(size_t i = 0;i < n * 60;i++){
				num *= 10;
				chars.push_back('0' + (num % 10));
			}
		}
		return std::string(chars.begin(), chars.end());
		
	}
	friend std::ostream& operator<<(std::ostream& out, const fixed<m, n, int_t>& o){
		for(size_t i = 0;i < m + n;i++){
			out << o.bits[i];
			if(i < m + n - 1){
				out << ", ";
			}
		}
		return out;
	}
	const int_t& operator[](size_t i)const{
		return bits[i];
	}
	int_t& operator[](size_t i){
		return bits[i];
	}
	int_t lzcnt(){
		int_t cnt = 0;
		size_t i = 0;
		while(bits[i] == 0)i++;
		return i*sizeof(int_t) * CHAR_BIT + lzcnt_intrinsic(bits[i]);
	}
	int_t tzcnt(){
		int_t cnt = 0;
		size_t i = m + n;
		while(bits[--i] == 0);
		return (m + n - 1 - i)*sizeof(int_t) * CHAR_BIT + tzcnt_intrinsic(bits[i]);
	}
	int_t popcnt(){
		int_t acc = 0;
		for(size_t i = 0;i < m+n;i++){
			acc += popcnt_intrinsic(bits[i]);
		}
		return acc;
	}
	bool operator==(const fixed<m, n>& o)const{
		for(size_t i = 0;i < m + n;i++){
			if(bits[i] != o[i])return false;
		}
		return true;
	}
	bool operator!=(const fixed<m, n>& o)const{
		for(size_t i = 0;i < m + n;i++){
			if(bits[i] != o[i])return true;
		}
		return false;
	}
	bool operator>(const fixed<m, n>& o)const{
		for(size_t i = 0;i < m + n - 1;i++){
			if(o[i] > bits[i])return false;
			if(bits[i] > o[i])return true;
		}
		return bits.back() > o.bits.back();
	}
	bool operator>=(const fixed<m, n>& o)const{
		for(size_t i = 0;i < m + n - 1;i++){
			if(o[i] > bits[i])return false;
			if(bits[i] > o[i])return true;
		}
		return bits.back() >= o.bits.back();
	}
	bool operator<(const fixed<m, n>& o)const{
		for(size_t i = 0;i < m + n - 1;i++){
			if(o[i] < bits[i])return false;
			if(bits[i] < o[i])return true;
		}
		return bits.back() < o.bits.back();
	}
	bool operator<=(const fixed<m, n>& o)const{
		for(size_t i = 0;i < m + n - 1;i++){
			if(o[i] < bits[i])return false;
			if(bits[i] < o[i])return true;
		}
		return bits.back() <= o.bits.back();
	}
	fixed operator<<(int64_t c)const{
		fixed a(*this);
		a.shiftLeft(c);
		return a;
	}
	fixed operator>>(int64_t c)const{
		fixed a(*this);
		a.shiftRight(c);
		return a;
	}
	fixed& operator<<=(int64_t c){
		shiftLeft(c);
		return *this;
	}
	fixed& operator>>=(int64_t c){
		shiftRight(c);
		return *this;
	}
	void shiftLeft(const int64_t c){
		if(c < 0){shiftRight(-c);return;}
		if(c == 0)return;
		if(c >= (m + n) * sizeof(int_t) * CHAR_BIT){
			std::fill(begin(), end(), 0);
			return;
		}
		const int64_t cs = c / (CHAR_BIT * sizeof(int_t));
		const int64_t bs = c % (CHAR_BIT * sizeof(int_t));
		for(size_t i = 0;i < m + n - cs - 1;i++){
			bits[i] = bits[i + cs] << bs | bits[i + cs + 1] >> (CHAR_BIT * sizeof(int_t) - bs);
		}
		bits[m + n - cs - 1] = bits.back() << bs;
		for(size_t i = m + n - cs; i < m + n;i++){
			bits[i] = 0;
		}
	}
	void shiftRight(const int64_t c){
		if(c < 0){shiftLeft(-c);return;}
		if(c == 0)return;
		if(c >= (m + n) * sizeof(int_t) * CHAR_BIT){
			std::fill(begin(), end(), 0);
			return;
		}
		const int64_t cs = c / (CHAR_BIT * sizeof(int_t));
		const int64_t bs = c % (CHAR_BIT * sizeof(int_t));
		for(int64_t i = m + n - 1;i >= cs + 1;i--){
			bits[i] = bits[i - cs] >> bs | bits[i - cs - 1] << (CHAR_BIT * sizeof(int_t) - bs);
		}
		bits[cs] = bits.front() >> bs;
		for(size_t i = 0; i < cs;i++){
			bits[i] = 0;
		}
	}
	constexpr auto begin(){return bits.begin();}
	constexpr auto cbegin()const{return bits.cbegin();}
	constexpr auto rbegin(){return bits.rbegin();}
	constexpr auto crbegin()const{return bits.crbegin();}
	constexpr auto   end(){return bits.  end();}
	constexpr auto  cend()const{return bits. cend();}
	constexpr auto  rend(){return bits. rend();}
	constexpr auto crend()const{return bits.crend();}
	constexpr auto begin()const{return bits.begin();}
	constexpr auto rbegin()const{return bits.rbegin();}
	constexpr auto   end()const{return bits.  end();}
	constexpr auto  rend()const{return bits. rend();}
	std::string bitstring(){
		std::string ret(bits.size() * sizeof(int_t) * CHAR_BIT + 1,'0');
		size_t adder = 0;
		for(size_t i = 0;i < bits.size();i++){
			for(size_t j = 1;j <= CHAR_BIT * sizeof(int_t);j++){
				ret[i * CHAR_BIT * sizeof(int_t) + j - 1 + adder] = '0' + !!(bits[i] & (int_t(1) << (CHAR_BIT * sizeof(int_t) - j)));
			}
			if(i == m - 1){
				adder = 1;
				ret[(i + 1) * CHAR_BIT * sizeof(int_t)] = '.';
			}
		}
		return ret;
	}
};
#endif
