#include "huffman.hpp"
#include "exceptions.hpp"
#include <queue>

Node::Node(int data, unsigned freq) : data(data), freq(freq), left(nullptr), right(nullptr) {}

std::map<int, std::string> HuffmanTree::encode(const std::vector<int>& data, const std::vector<int>& freq) {
    struct Compare {
        bool operator()(const std::unique_ptr<Node>& l, const std::unique_ptr<Node>& r) const {
            return l->freq > r->freq;
        }
    };

    std::priority_queue<std::unique_ptr<Node>, std::vector<std::unique_ptr<Node>>, Compare> minHeap;

    // Create leaf nodes
    for (size_t i = 0; i < data.size(); ++i) {
        if (freq[i] > 0) {
            minHeap.push(std::make_unique<Node>(data[i], freq[i]));
        }
    }

    // Build Huffman tree
    while (minHeap.size() > 1) {
        auto left = std::move(const_cast<std::unique_ptr<Node>&>(minHeap.top()));
        minHeap.pop();
        auto right = std::move(const_cast<std::unique_ptr<Node>&>(minHeap.top()));
        minHeap.pop();

        auto parent = std::make_unique<Node>(-1, left->freq + right->freq);
        parent->left = std::move(left);
        parent->right = std::move(right);
        minHeap.push(std::move(parent));
    }

    root = std::move(const_cast<std::unique_ptr<Node>&>(minHeap.top()));
    std::map<int, std::string> huffmanCodes;
    generateCodes(root.get(), "", huffmanCodes);
    return huffmanCodes;
}

void HuffmanTree::generateCodes(Node* node, std::string code, std::map<int, std::string>& huffmanCodes) {
    if (!node) return;

    if (node->data != -1) {
        huffmanCodes[node->data] = code;
    }

    generateCodes(node->left.get(), code + "0", huffmanCodes);
    generateCodes(node->right.get(), code + "1", huffmanCodes);
}

int HuffmanTree::decode(const std::string& encodedStr, const std::map<int, std::string>& huffmanCodes) {
    Node* current = root.get();

    for (char bit : encodedStr) {
        if (!current) {
            throw ImageException("Invalid Huffman code sequence");
        }

        current = (bit == '0') ? current->left.get() : current->right.get();

        if (current && !current->left && !current->right) {
            return current->data;
        }
    }
    throw ImageException("Incomplete Huffman code sequence");
}
