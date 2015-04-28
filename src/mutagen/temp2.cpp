class Thing
{
public:
    Thing(const int n)  { p = new long[n]; }
    long* p;
};

void foo()
{
    Thing b(16);
    const long * bp = b.p;
    const long r = bp[0];
}
