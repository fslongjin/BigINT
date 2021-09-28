//
// Created by longjin on 2021/9/25.
//

//#pragma GCC optimize("O2")

#include "Bigint.h"
#include<iostream>
#include<stdexcept>
#include "stack"

#define LARGER 1
#define EQUAL 0
#define SMALLER -1

using namespace std;

Bigint::Bigint() {
    this->begin = new node();
    this->end = this->begin;
    this->begin->num = 0;
}


Bigint::Bigint(long long x) {
    this->begin = new node();
    this->end = this->begin;

    //防止由于x为0而造成异常
    if (x == 0)
        this->begin->num = 0;
    else
        this->translate_num2BigInt(x);

}

Bigint::Bigint(std::string x) {
    this->begin = new node();
    this->end = this->begin;

    //防止由于x为0而造成异常
    //this->begin->num = 0;
    this->translate_string2BigInt(x);
}

void Bigint::translate_num2BigInt(long long x) {
    //把数字转为bigInt

    if (x < 0) {
        this->sign = false;
        x = -x;
    } else this->sign = true;

    while (x > 0) {
        if (this->begin == this->end && this->begin->num == -1)//该位未初始化
            this->begin->num = x % 10;
        else {
            new_node_link_to_begin();
            this->begin->num = x % 10;
        }
        x /= 10;
    }

}

//生成新节点并连接到begin前面
void Bigint::new_node_link_to_begin() {
    this->begin->prev = new node();
    this->begin = this->begin->prev;
    ++size;
}


//重载输出
std::ostream &operator<<(ostream &out, const Bigint &p) {
    node *nd = p.end;
    string ans;
    while (nd != nullptr) {
        ans = to_string(nd->num) + ans;
        nd = nd->prev;
    }
    if (!p.sign)
        ans = "-" + ans;
    out << ans;

    return out;
}

void clear_num(Bigint *p) {
//删除所有节点

    while (p->begin != p->end) {
        node *tmp = p->end;
        p->end = p->end->prev;
        free(tmp);
    }


    p->size = 1;
    p->sign = true;
    p->end->num = 0;
}

Bigint &Bigint::operator=(const Bigint &p) {


    while (this->begin != this->end) {
        node *tmp = this->end;
        this->end = this->end->prev;
        free(tmp);
    }

    size = 1;
    this->sign = p.sign;
    this->end->num = p.end->num;
    node *ptr = p.end->prev;
    while (ptr != nullptr) {
        new_node_link_to_begin();
        this->begin->num = ptr->num;
        ptr = ptr->prev;

    }
    return *this;

}

Bigint::~Bigint() {
    node *ptr = this->end;
    node *tmp;
    while (ptr != nullptr) {
        tmp = ptr;
        ptr = ptr->prev;
        free(tmp);
    }
}

std::istream &operator>>(istream &in, Bigint &p) {
    clear_num(&p);
    string s;
    in >> s;
    p.translate_string2BigInt(s);
    return in;
}

void Bigint::translate_string2BigInt(std::string x) {
    this->begin->num = x[x.length() - 1] - '0';
    int end_pos = 0;
    if (x[0] == '-') {
        this->sign = false;
        end_pos = 1;
    } else if (x[0] == '+') {
        this->sign = true;
        end_pos = 1;
    } else if (x[0] - '0' > 9 || x[0] < '0') {
        throw invalid_argument("syntax error");
    }

    this->begin->num = x[x.length() - 1] - '0';

    for (long long i = x.length() - 2; i >= end_pos; --i) {
        if (x[i] - '0' > 9 || x[i] < '0') {
            throw invalid_argument("syntax error");
        }
        new_node_link_to_begin();
        this->begin->num = x[i] - '0';
    }
}

//拷贝构造函数
Bigint::Bigint(const Bigint &p) {
    this->begin = new node();
    this->end = this->begin;
    this->begin->num = 0;

    this->sign = p.sign;
    this->size = p.size;
    this->begin->num = p.end->num;

    node *ndp = p.end->prev;

    while (ndp != nullptr) {
        new_node_link_to_begin();
        this->begin->num = ndp->num;
        ndp = ndp->prev;
    }

}

/**
  * 清除前导零，更新size
*/
void Bigint::clean_leading_zero() {
    node *nd = this->end;

    stack<node *> prev;
    while (nd != this->begin) {
        prev.push(nd);
        nd = nd->prev;
    }
    while (!prev.empty() && size > 1) {
        if (nd->num == 0) {
            this->begin = prev.top();
            this->begin->prev = nullptr;
            prev.pop();
            --size;
            free(nd);
            nd = this->begin;
        } else break;
    }
    //防止 -0 的情况出现
    if (size == 1)
        if (this->end->num == 0)
            this->sign = true;

}

Bigint Bigint::operator+(const Bigint &p) const {

    Bigint large;
    Bigint small;
    if (*this < p) {
        large = p;
        small = *this;
    } else {
        large = *this;
        small = p;
    }

    if (this->sign == p.sign) {
        //符号相同，逐位相加
        node *nd = large.end;
        node *ndp = small.end;
        int carry = 0;
        while (nd != nullptr && ndp != nullptr) {
            nd->num = nd->num + ndp->num + carry;
            if (nd->num >= 10)
                carry = nd->num / 10;
            else carry = 0;

            nd->num %= 10;

            nd = nd->prev;
            ndp = ndp->prev;
        }

        //将突出来的数位加上
        while (nd == nullptr && ndp != nullptr) {
            large.new_node_link_to_begin();
            nd = large.begin;
            nd->num = ndp->num + carry;
            if (nd->num >= 10)
                carry = nd->num / 10;
            else carry = 0;

            nd->num %= 10;

            nd = nd->prev;
            ndp = ndp->prev;

        }


        if (carry) {
            //进位不为0，需要继续进位
            large.new_node_link_to_begin();
            large.begin->num = carry;
        }

        large.clean_leading_zero();
        return Bigint(large);

    } else {
        //符号不相同

        small.sign = true;
        Bigint ans = large - small;
        ans.clean_leading_zero();
        return ans;
    }

}

Bigint Bigint::operator-(const Bigint &p) const {


    int borrow = 0;
    Bigint large;
    Bigint small;
    large = *this;
    small = p;


    //只计算两个正数之间的相减
    if (large.sign == small.sign && large.sign) {
        if (large < small) {
            //被减数比减数小
            large = p;
            small = *this;
            large.sign = false;
        }
        node *nd = large.end;
        node *ndp = small.end;

        while (nd != nullptr && ndp != nullptr) {
            nd->num = nd->num - ndp->num - borrow;
            if (nd->num < 0) {
                borrow = 1;
                nd->num += 10;
            } else borrow = 0;

            nd = nd->prev;
            ndp = ndp->prev;
        }
        while (borrow && nd != nullptr) {
            nd->num -= borrow;
            if (nd->num < 0) {
                borrow = 1;
                nd->num += 10;
            } else borrow = 0;
            nd = nd->prev;
        }
        large.clean_leading_zero();
        return Bigint(large);

    } else if (large.sign == small.sign && !large.sign) {
        //两个负数相减
        small.sign = true;
        return large + small;
    } else {
        //一正一负
        if (!small.sign) {
            //正减负
            small.sign = true;
            return large + small;
        } else {
            //负减正
            large.sign = true;
            Bigint ans = large + small;
            ans.sign = false;
            return ans;
        }
    }


}

bool Bigint::operator<(const Bigint &p) const {

    int ans = compare_with((*this), p);
    if (ans == SMALLER)
        return true;
    else return false;
    /*
    if (this->sign < p.sign)
        return true;

    else if (this->sign > p.sign)
        return false;

    //符号相同
    if (this->size < p.size) {
        if (!this->sign)
            return false;
        else return true;
    } else if (this->size > p.size) {
        if (!this->sign)
            return true;
        else return false;

    }

    node *nd = this->end;
    node *ndp = p.end;
    while (nd != nullptr && ndp != nullptr) {
        if (nd->num > ndp->num)
            return false;
        else if (nd->num < ndp->num)
            return true;

        nd = nd->prev;
        ndp = ndp->prev;
    }
    //相等
    if (nd == nullptr && ndp == nullptr)
        return false;
    else if (nd == nullptr && ndp != nullptr)
        return true;
        //大于
    else if (nd != nullptr && ndp == nullptr)
        return false;
    return false;
*/
}

Bigint Bigint::operator*(const Bigint &p) const {
    Bigint ans;
    if (this->sign == p.sign)
        ans.sign = true;
    else ans.sign = false;

    node *nd = this->end;

    int x;
    long long i = 0;

    while (nd != nullptr) {
        x = 0;
        node *ndp = p.end;
        long long j = 0;

        //获取到ans的第i位的指针
        node *nda = ans.ith_pointer(i);

        while (ndp != nullptr) {
            nda->num += nd->num * ndp->num + x;
            x = nda->num / 10;
            nda->num %= 10;

            ndp = ndp->prev;

            if (nda->prev == nullptr) {
                ans.new_node_link_to_begin();
                nda->prev->num = 0;
            }
            nda = nda->prev;
            ++j;
        }
        nda->num = x;

        nd = nd->prev;
        ++i;
    }
    ans.clean_leading_zero();
    return Bigint(ans);
}

/**
    * 获取第i位的指针
    * @param i 位数
    */
node *Bigint::ith_pointer(long long int i) {
    node *nda = this->end;
    for (int k = 0; k <= i; ++k) {
        if (nda->prev == nullptr) {
            this->new_node_link_to_begin();
            this->begin->num = 0;
        }
        if (k != i)
            nda = nda->prev;
    }
    return nda;
}


int Bigint::compare_with(const Bigint &a, const Bigint &p) {
    if (a.sign < p.sign)
        return SMALLER;

    else if (a.sign > p.sign)
        return LARGER;

    //符号相同
    if (a.size < p.size) {
        if (!a.sign)
            return LARGER;
        else return SMALLER;
    } else if (a.size > p.size) {
        if (!a.sign)
            return SMALLER;
        else return LARGER;

    }

    //位数相同，符号相同

    node *nd = a.end;
    node *ndp = p.end;

    stack<node *> prev_nd;
    stack<node *> prev_ndp;

    while (nd != nullptr) {
        prev_nd.push(nd);
        nd = nd->prev;
    }
    while (ndp != nullptr) {
        prev_ndp.push(ndp);
        ndp = ndp->prev;
    }

    while (!prev_nd.empty()) {
        nd = prev_nd.top();
        prev_nd.pop();

        ndp = prev_ndp.top();
        prev_ndp.pop();

        if (a.sign) {
            //正数
            if (nd->num > ndp->num)
                return LARGER;
            else if (nd->num < ndp->num)
                return SMALLER;
        } else {
            //负数
            if (nd->num > ndp->num)
                return SMALLER;
            else if (nd->num < ndp->num)
                return LARGER;
        }
    }

    //相等
    if (nd == a.end && ndp == p.end)
        return EQUAL;
    else {
        //不应该到达这里
        throw invalid_argument("At Compare: final pointer is not a nullptr");
    }
}

bool Bigint::operator==(const Bigint &p) const {
    int ans = compare_with(*this, p);
    if (ans == EQUAL)
        return true;
    else return false;
}


Bigint Bigint::operator/(const Bigint &p) const {
    if (this->sign == p.sign) {
        if (this->sign && (*this) < p)
            return Bigint();
        else if (!this->sign && p < (*this))
            return Bigint();
    } else if ((*this) == 0)
        return Bigint();
    else if ((*this) == p)
        return Bigint(1);
    else if (p == 0)
        //除零错误
        throw invalid_argument("Div 0 error!");

    string a = "", b = "";

    node *nda = this->end;

    while (nda != nullptr) {
        a = to_string(nda->num) + a;
        nda = nda->prev;
    }

    node *ndb = p.end;

    while (ndb != nullptr) {
        b = to_string(ndb->num) + b;
        ndb = ndb->prev;
    }

    int l_pos = 0, r_pos = b.length() - 1;

    Bigint ans = Bigint();
    if (this->sign == p.sign)
        ans.sign = true;
    else ans.sign = false;


    while (true) {
        if (r_pos >= a.length())
            break;
        if (!check_div(a, b, l_pos, r_pos)) {
            ++r_pos;
            continue;
        } else {
            int tmp_ans = 0;
            while (true) {
                if (!check_div(a, b, l_pos, r_pos)) {
                    break;
                }

                //能够减
                int js = b.length() - 1;
                for (int i = r_pos; i >= l_pos && js >= 0; --i) {
                    int tmp = (a[i] - '0') - (b[js] - '0');

                    if (tmp < 0) {
                        a[i - 1] -= 1;

                        a[i] = '0' + tmp + 10;
                    } else {
                        a[i] = '0' + tmp;
                    }
                    --js;
                    if ((a[l_pos] - '0') == 0)
                        ++l_pos;
                }
                ++tmp_ans;

            }
            int js = 1;
            int to = a.length() - r_pos;
            node *nd = ans.end;
            while (js < to) {
                if (nd->prev == nullptr) {
                    ans.new_node_link_to_begin();
                    ans.begin->num = 0;
                }
                nd = nd->prev;
                ++js;
            }
            nd->num = tmp_ans;

            //更新l_pos
            for (int i = l_pos; i <= r_pos; ++i) {
                if (a[i] == '0')
                    ++l_pos;
            }

        }


    }
    ans.clean_leading_zero();
    return ans;
}

bool Bigint::check_div(const string &a, const string &b, const int &l_pos, const int &r_pos) {
    long long d = r_pos - l_pos + 1;
    if (d > b.length())
        return true;
    else if (d < b.length())
        return false;

    for (int i = 0; i < b.length(); ++i) {
        if (a[l_pos + i] < b[i])
            return false;
        else if (a[l_pos + i] > b[i])
            return true;
    }
    return true;
}

bool Bigint::operator<=(const Bigint &p) const {
    int ans = compare_with((*this), p);
    if (ans == SMALLER || ans == EQUAL)
        return true;
    else return false;
}

bool Bigint::operator>(const Bigint &p) const {
    int ans = compare_with((*this), p);
    if (ans == LARGER)
        return true;
    else return false;
}

bool Bigint::operator>=(const Bigint &p) const {
    int ans = compare_with((*this), p);
    if (ans == LARGER || ans == EQUAL)
        return true;
    else return false;
}

Bigint pow(const Bigint &x, const Bigint n) {
    Bigint ans = Bigint(x);

    for (Bigint i = Bigint(1); i < n; i = i + 1) {
        ans = ans * x;
    }
    return ans;
}

Bigint qpow(const Bigint &x, long long n) {
    Bigint ans = Bigint(1);
    Bigint xx = Bigint(x);

    while (n) {
        if (n & 1)
            ans = ans * xx;
        xx = xx * xx;
        n >>= 1;
    }
    return ans;

}