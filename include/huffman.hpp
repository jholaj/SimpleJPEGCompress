#pragma once

#include <memory>
#include <map>
#include <string>
#include <vector>
#include "pixel_block.hpp"

struct Node {
    int data;
    unsigned freq;
    std::unique_ptr<Node> left, right;

    Node(int data, unsigned freq);
};

class HuffmanTree {
public:
    std::map<int, std::string> encode(const std::vector<int>& data, const std::vector<int>& freq);
    int decode(const std::string& encodedStr, const std::map<int, std::string>& huffmanCodes);

private:
    std::unique_ptr<Node> root;
    void generateCodes(Node* node, std::string code, std::map<int, std::string>& huffmanCodes);
};
