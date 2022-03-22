#pragma once
#include "structures.h"

My_vector operator+ (const My_vector& lhs, const My_vector& rhs) {
	return { lhs.GetX() + rhs.GetX(), lhs.GetY() + rhs.GetY(), lhs.GetZ() + rhs.GetZ() };
}
My_vector operator- (const My_vector& lhs, const My_vector& rhs) {
	return { lhs.GetX() - rhs.GetX(), lhs.GetY() - rhs.GetY(), lhs.GetZ() - rhs.GetZ() };
}
My_vector operator* (const My_vector& vec, const double a) {
	return { vec.GetX() * a, vec.GetY() * a, vec.GetZ() * a };
}
My_vector operator* (const double a, const My_vector& vec) {
	return { vec.GetX() * a, vec.GetY() * a, vec.GetZ() * a };
}
My_vector operator/ (const My_vector& vec, const double a) {
	return { vec.GetX() / a, vec.GetY() / a, vec.GetZ() / a };
}
My_vector operator/ (const My_vector& lhs, const My_vector& rhs) {
	return { lhs.GetY() * rhs.GetZ() - lhs.GetZ() * rhs.GetY(), lhs.GetZ() * rhs.GetX() - lhs.GetX() * rhs.GetZ(), lhs.GetX() * rhs.GetY() - lhs.GetY() * rhs.GetX() };
}
std::pair<My_vector, My_vector> operator* (const std::pair<My_vector, My_vector> vec, const double a) {
	My_vector f = vec.first * a;
	My_vector s = vec.second * a;
	return { f, s };
}
std::pair<My_vector, My_vector> operator* (const double a, const std::pair<My_vector, My_vector> vec) {
	My_vector f = vec.first * a;
	My_vector s = vec.second * a;
	return { f, s };
}
std::pair<My_vector, My_vector> operator+ (const std::pair<My_vector, My_vector> lhs, const std::pair<My_vector, My_vector> rhs) {
	My_vector f = lhs.first + rhs.first;
	My_vector s = lhs.second + rhs.second;
	return { f, s };
}

std::vector<My_vector> operator* (const std::vector<My_vector> vec, const double a) {
	std::vector<My_vector> tmp;
	for (const auto& x : vec) {
		tmp.push_back(x * a);
	}
	return tmp;
}
std::vector<My_vector> operator* (const double a, const std::vector<My_vector> vec) {
	std::vector<My_vector> tmp;
	for (const auto& x : vec) {
		tmp.push_back(x * a);
	}
	return tmp;
}
std::vector<My_vector> operator+ (const std::vector<My_vector> lhs, const std::vector<My_vector> rhs) {
	std::vector<My_vector> tmp;
	if (lhs.size() == rhs.size()) {
		for (size_t i = 0; i < lhs.size(); i++) {
			tmp.push_back(lhs[i] + rhs[i]);
		}
	}
	else throw std::out_of_range("sizes of vectors are not equal");
	return tmp;
}
