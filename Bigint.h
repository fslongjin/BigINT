//
// Created by longjin on 2021/9/25.
//
//#pragma GCC optimize("O2")

#ifndef INC_BIGINT_BIGINT_H
#define INC_BIGINT_BIGINT_H


#include <string>
#include <unordered_map>
#include "iostream"

struct node {
    int num;
    node *prev;

    node() {
        prev = nullptr;
        num = -1;
    }
};


class Bigint {
public:
    Bigint();


    Bigint(long long);

    Bigint(std::string);

    Bigint(Bigint const &);

    friend std::ostream &operator<<(std::ostream &out, const Bigint &p);

    friend std::istream &operator>>(std::istream &in, Bigint &p);


    //使得Bigint清零
    friend void clear_num(Bigint *);

    Bigint &operator=(const Bigint &p);


    Bigint operator+(const Bigint &p) const;

    Bigint operator-(const Bigint &p) const;
    Bigint operator*(const Bigint &p)const;
    Bigint operator/(const Bigint &p)const;




    bool operator<(const Bigint &p) const;
    bool operator<=(const Bigint &p) const;
    bool operator>(const Bigint &p) const;
    bool operator>=(const Bigint &p) const;
    bool operator==(const Bigint &p)const;

    ~Bigint();

protected:
    void new_node_link_to_begin();

private:
    node *begin;
    node *end;
    int size = 1;
    bool sign = true; //正为1，负为0
    void translate_num2BigInt(long long);

    void translate_string2BigInt(std::string);

    /**
     * 清除前导零，更新size
     */
    void clean_leading_zero();

    /**
     * 获取第i位的指针
     * @param i 位数
     */
    node* ith_pointer(long long i);

    /**
     * 比较a和p的大小
     * @param a
     * @param p
     * @return
     */
    static int compare_with(const Bigint& a, const Bigint& p);

    //检查能否在此次除法中相减
    static bool check_div(const std::string &a, const std::string &b, const int &l_pos, const int& r_pos);


};

//乘方
Bigint pow(const Bigint& x, const Bigint n);
//快速幂

Bigint qpow(const Bigint &x, long long n);


#endif //INC_BIGINT_BIGINT_H
