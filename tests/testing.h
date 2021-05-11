#ifndef TESTING_H
#define TESTING_H

#include <map>
#include <string>
#include <vector>

#include "dist.h"
#include "FragmentModels/FragmentIon.h"

void serializeDist(std::map<std::string, std::vector<float> > dist, std::string filename);

std::vector<FragmentIon> flattenFragments(const std::vector< std::vector<FragmentIon> > &in);

void writeFragments(const std::vector<FragmentIon> &in, std::string filename);

std::string serializeFragments(const std::vector<FragmentIon> &in);

std::vector<FragmentIon> unserializeFragments(std::string filename);

bool compareFragments(const std::vector<FragmentIon> &a, const std::vector<FragmentIon> &b);

void printDifference(std::vector<FragmentIon> a, std::vector<FragmentIon> b);

bool isEqual(double x, double y, int precision);

#endif // TESTING_H
