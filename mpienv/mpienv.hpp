namespace MPIEnv
{
    int initialize(int *argc, char **argv[])
    {
        if (is_initialized()) return MPI_ERR_OTHER;
        cache.reset();
        return MPI_Init(argc, argv);
    }

    int finalize()
    {
        if (is_finalized()) return MPI_SUCCESS;
        cache.reset();
        return MPI_Finalize();
    }

    void exit(int err)
    {
        finalize();
        std::exit(err);
    }

    bool is_initialized()
    {
        int flag;
        MPI_Initialized(&flag);
        return (flag != 0);
    }

    bool is_finalized()
    {
        int flag;
        MPI_Finalized(&flag);
        return (flag != 0);
    }

    bool comms_equal(MPI_Comm lhs, MPI_Comm rhs)
    {
        int result;
        MPI_Comm_compare(lhs, rhs, &result);
        return (result == MPI_IDENT) || (result == MPI_CONGRUENT);
    }

    bool type_info_compare::operator()(const std::type_info *lhs, const std::type_info *rhs) const
    {
        return lhs->before(*rhs);
    }

    template <class T> void mpi_create_type(MPI_Datatype& dtype)
    {
        if constexpr (is_array_type<T>)
        {
            using U = typename array_info<T>::value_type;
            static constexpr int count = array_info<T>::size;
            MPI_Type_contiguous(count, mpi_type<U>(), &dtype);
        }
        else
        {
            MPI_Type_contiguous(sizeof(T), MPI_CHAR, &dtype);
        }

        MPI_Type_commit(&dtype);
    }

    template <class T> MPI_Datatype mpi_commit_type()
    {
        const std::type_info *t = &typeid(T);
        MPI_Datatype dtype = cache.get_type(t);

        if (dtype == MPI_DATATYPE_NULL)
        {
            mpi_create_type<T>(dtype);
            cache.set_type(t, dtype);
        }

        return dtype;
    }

    template <class T> MPI_Datatype mpi_type()
    {
        if      constexpr (std::same_as<T, char>)               return MPI_CHAR;
        else if constexpr (std::same_as<T, signed char>)        return MPI_SIGNED_CHAR;
        else if constexpr (std::same_as<T, short>)              return MPI_SHORT;
        else if constexpr (std::same_as<T, int>)                return MPI_INT;
        else if constexpr (std::same_as<T, long>)               return MPI_LONG;
        else if constexpr (std::same_as<T, long long>)          return MPI_LONG_LONG;
        else if constexpr (std::same_as<T, unsigned char>)      return MPI_UNSIGNED_CHAR;
        else if constexpr (std::same_as<T, unsigned short>)     return MPI_UNSIGNED_SHORT;
        else if constexpr (std::same_as<T, unsigned int>)       return MPI_UNSIGNED;
        else if constexpr (std::same_as<T, unsigned long>)      return MPI_UNSIGNED_LONG;
        else if constexpr (std::same_as<T, unsigned long long>) return MPI_UNSIGNED_LONG_LONG;
        else if constexpr (std::same_as<T, float>)              return MPI_FLOAT;
        else if constexpr (std::same_as<T, double>)             return MPI_DOUBLE;
        else if constexpr (std::same_as<T, long double>)        return MPI_LONG_DOUBLE;
        else if constexpr (std::same_as<T, bool>)               return MPI_CXX_BOOL;
        else                                                    return mpi_commit_type<T>();
    }

    TypeCache::TypeCache() {}
    TypeCache::~TypeCache() { reset(); }

    MPI_Datatype TypeCache::get_type(const std::type_info* t)
    {
        auto pos = type_map.find(t);

        if (pos == type_map.end())
            return MPI_DATATYPE_NULL;
        else
            return pos->second;
    }

    void TypeCache::set_type(const std::type_info *t, MPI_Datatype dtype)
    {
        type_map.insert({t, dtype});
    }

    void TypeCache::free_types()
    {
        if (!is_finalized())
            for (auto& [t, dtype] : type_map)
                MPI_Type_free(&dtype);
    }

    void TypeCache::reset() { free_types(); }

    void Comm::init(MPI_Comm comm)
    {
        MPI_Comm_dup(comm, commbuf);
        MPI_Comm_rank(commbuf[0], &myrank);
        MPI_Comm_size(commbuf[0], &nprocs);
    }

    Comm::Comm() : Comm(MPI_COMM_NULL) {}

    Comm::Comm(MPI_Comm comm) { init(comm); }

    Comm::Comm(const Comm& rhs) { init(rhs.getcomm()); }

    Comm::~Comm() { MPI_Comm_free(commbuf); }

    void Comm::swap(Comm& rhs) noexcept
    {
        MPI_Comm tmp[1];

        std::swap(myrank, rhs.myrank);
        std::swap(nprocs, rhs.nprocs);

        std::memcpy(tmp, commbuf, sizeof(MPI_Comm));
        std::memcpy(commbuf, rhs.commbuf, sizeof(MPI_Comm));
        std::memcpy(rhs.commbuf, tmp, sizeof(MPI_Comm));
    }

    bool Comm::operator==(const Comm &rhs) const
    {
        return comms_equal(getcomm(), rhs.getcomm());
    }

    Comm& Comm::operator=(const Comm& rhs)
    {
        init(rhs.getcomm());
        return *this;
    }

    int Comm::barrier() const
    {
        return MPI_Barrier(getcomm());
    }

    template <class T>
    bool Comm::is_same_val(T val) const
    {
        auto comm = getcomm();
        bool is_same;

        if (nprocs == 1) return true;

        auto dtype = mpi_type<T>();

        if constexpr (std::is_arithmetic_v<T>)
        {
            T buf[2] = {-val, val};
            MPI_Allreduce(MPI_IN_PLACE, buf, 2, dtype, MPI_MIN, comm);
            is_same = (buf[0] == -buf[1]);
        }
        else
        {
            std::vector<T> buffer;

            if (myrank == 0) buffer.resize(nprocs);

            MPI_Gather(&val, 1, dtype, buffer.data(), 1, dtype, 0, comm);

            if (myrank == 0) is_same = std::equal(buffer.begin()+1, buffer.end(), buffer.begin());

            MPI_Bcast(&is_same, 1, mpi_type<bool>(), 0, comm);
        }

        return is_same;
    }

    template <class T>
    bool Comm::are_same_vals(const std::vector<T>& vals) const
    {
        auto comm = getcomm();

        if (nprocs == 1) return true;

        int count = vals.size();

        if (!is_same_val(count)) return false;

        if (count == 0) return true;

        auto dtype = mpi_type<T>();

        std::vector<T> buf = vals, rbuf;

        if constexpr (std::is_arithmetic_v<T>)
        {
            buf.resize(count<<1);
            std::transform(buf.begin(), buf.begin() + count, buf.begin() + count, [](T x) { return -x; });

            MPI_Allreduce(MPI_IN_PLACE, buf.data(), count<<1, dtype, MPI_MIN, comm);

            for (int i = 0; i < count; ++i)
                if (-buf[i] != buf[i+count])
                    return false;

            return true;
        }
        else
        {
            if (myrank == 0) rbuf.resize(count*nprocs);

            MPI_Gather(buf.data(), count, dtype, rbuf.data(), count, dtype, 0, comm);

            bool are_same = true;

            if (myrank == 0)
            {
                auto it = rbuf.begin() + count;

                for (int i = 1; are_same && i < nprocs; it += count, ++i)
                    if (!std::equal(buf.begin(), buf.end(), it))
                        are_same = false;
            }

            MPI_Bcast(&are_same, 1, mpi_type<bool>(), 0, comm);
            return are_same;
        }
    }

    template <class T>
    int Comm::reduce(const T* sendbuf, T* recvbuf, int count, int root, MPI_Op op) const
    {
        return MPI_Reduce(sendbuf, recvbuf, count, mpi_type<T>(), op, root, getcomm());
    }

    template <class T>
    int Comm::reduce(const T& sendbuf, T& recvbuf, int root, MPI_Op op) const
    {
        return reduce(&sendbuf, &recvbuf, 1, root, op);
    }

    template <class T>
    int Comm::reduce(const std::vector<T>& sendbuf, std::vector<T>& recvbuf, int root, MPI_Op op) const
    {
        int count = sendbuf.size();
        if (!is_same_val(count)) return MPI_ERR_OTHER;
        if (myrank == root) recvbuf.resize(count);
        return reduce(sendbuf.data(), recvbuf.data(), count, root, op);
    }

    template <class T>
    int Comm::bcast(T* buffer, int count, int root) const
    {
        return MPI_Bcast(buffer, count, mpi_type<T>(), root, getcomm());
    }

    template <class T>
    int Comm::bcast(T& buffer, int root) const
    {
        return bcast(&buffer, 1, root);
    }

    template <class T>
    int Comm::bcast(std::vector<T>& buffer, int root) const
    {
        int count = buffer.size();
        MPI_Bcast(&count, 1, MPI_INT, root, getcomm());
        buffer.resize(count);
        return bcast(buffer.data(), count, root);
    }

    template <class T>
    int Comm::exscan(const T* sendbuf, T* recvbuf, int count, MPI_Op op, T identity) const
    {
        std::fill(recvbuf, recvbuf + count, identity);
        return MPI_Exscan(sendbuf, recvbuf, count, mpi_type<T>(), op, getcomm());
    }

    template <class T>
    int Comm::exscan(const T& sendbuf, T& recvbuf, MPI_Op op, T identity) const
    {
        return exscan(&sendbuf, &recvbuf, 1, op, identity);
    }

    template <class T> int Comm::allreduce(const T* sendbuf, T* recvbuf, int count, MPI_Op op) const
    {
        return MPI_Allreduce(sendbuf, recvbuf, count, mpi_type<T>(), op, getcomm());
    }

    template <class T> int Comm::allreduce(const T& sendbuf, T& recvbuf, MPI_Op op) const
    {
        return allreduce(&sendbuf, &recvbuf, 1, op);
    }

    template <class T> int Comm::allreduce(T* buffer, int count, MPI_Op op) const
    {
        return MPI_Allreduce(MPI_IN_PLACE, buffer, count, mpi_type<T>(), op, getcomm());
    }

    template <class T> int Comm::allreduce(T& buffer, MPI_Op op) const
    {
        return allreduce(&buffer, 1, op);
    }

    template <class T> int Comm::allreduce(std::vector<T>& buffer, MPI_Op op) const
    {
        return allreduce(buffer.data(), static_cast<int>(buffer.size()), op);
    }

    template <class T> int Comm::allreduce(const std::vector<T>& sendbuf, std::vector<T>& recvbuf, MPI_Op op) const
    {
        int count = sendbuf.size();
        if (!is_same_val(count)) return MPI_ERR_OTHER;
        recvbuf.resize(count);
        return allreduce(sendbuf.data(), recvbuf.data(), count, op);
    }

    template <class T> int Comm::gather(const T* sendbuf, int count, T* recvbuf, int root) const
    {
        return MPI_Gather(sendbuf, count, mpi_type<T>(), recvbuf, count, mpi_type<T>(), root, getcomm());
    }

    template <class T> int Comm::gather(const T& sendbuf, std::vector<T>& recvbuf, int root) const
    {
        if (myrank == root) recvbuf.resize(nprocs);
        return gather(&sendbuf, 1, recvbuf.data(), root);
    }

    template <class T> int Comm::gather(const std::vector<T>& sendbuf, std::vector<T>& recvbuf, int root) const
    {
        int count = sendbuf.size();
        if (!is_same_val(count)) return MPI_ERR_OTHER;
        if (myrank == root) recvbuf.resize(count*nprocs);
        return gather(sendbuf.data(), count, recvbuf.data(), root);
    }

    template <class T> int Comm::allgather(const T* sendbuf, int count, T* recvbuf) const
    {
        return MPI_Allgather(sendbuf, count, mpi_type<T>(), recvbuf, count, mpi_type<T>(), getcomm());
    }

    template <class T> int Comm::allgather(T* buffer, int count) const
    {
        return MPI_Allgather(MPI_IN_PLACE, count, mpi_type<T>(), buffer, count, mpi_type<T>(), getcomm());
    }

    template <class T> int Comm::allgather(std::vector<T>& buffer) const
    {
        int count = buffer.size();
        if (!is_same_val(count)) return MPI_ERR_OTHER;
        if (count % nprocs != 0) return MPI_ERR_OTHER;
        return allgather(buffer.data(), count/nprocs);
    }

    template <class T> int Comm::gatherv(const std::vector<T>& sendbuf, std::vector<T>& recvbuf, int root) const
    {
        auto comm = getcomm();
        auto dtype = mpi_type<T>();

        int sendcount = sendbuf.size();
        std::vector<int> recvcounts, displs;

        if (myrank == root)
        {
            recvcounts.resize(nprocs);
            displs.resize(nprocs);
        }

        int err = MPI_Gather(&sendcount, 1, MPI_INT, recvcounts.data(), 1, MPI_INT, root, comm);
        if (err != MPI_SUCCESS) return err;

        if (myrank == root)
        {
            std::exclusive_scan(recvcounts.begin(), recvcounts.end(), displs.begin(), static_cast<int>(0));
            recvbuf.resize(recvcounts.back() + displs.back());
        }

        return MPI_Gatherv(sendbuf.data(), sendcount, dtype, recvbuf.data(), recvcounts.data(), displs.data(), dtype, root, comm);
    }

    template <class T> int Comm::allgatherv(const std::vector<T>& sendbuf, std::vector<T>& recvbuf) const
    {
        auto comm = getcomm();
        auto dtype = mpi_type<T>();

        std::vector<int> recvcounts(nprocs), displs(nprocs);

        recvcounts[myrank] = sendbuf.size();

        int err = MPI_Allgather(MPI_IN_PLACE, 1, MPI_INT, recvcounts.data(), 1, MPI_INT, comm);
        if (err != MPI_SUCCESS) return err;

        std::exclusive_scan(recvcounts.begin(), recvcounts.end(), displs.begin(), static_cast<int>(0));
        recvbuf.resize(recvcounts.back() + displs.back());

        return MPI_Allgatherv(sendbuf.data(), recvcounts[myrank], dtype, recvbuf.data(), recvcounts.data(), displs.data(), dtype, comm);
    }

    template <class T> int Comm::scatter(const T* sendbuf, int count, T* recvbuf, int root) const
    {
        return MPI_Scatter(sendbuf, count, mpi_type<T>(), recvbuf, count, mpi_type<T>(), root, getcomm());
    }

    template <class T> int Comm::scatter(const std::vector<T>& sendbuf, T& recvbuf, int root) const
    {
        int valid = !!(sendbuf.size() == nprocs);
        bcast(valid, root);
        if (!valid) return MPI_ERR_OTHER;
        return scatter(sendbuf.data(), 1, &recvbuf, root);
    }

    template <class T> int Comm::scatter(const std::vector<T>& sendbuf, std::vector<T>& recvbuf, int root) const
    {
        int count = sendbuf.size();
        int valid = !!(count % nprocs == 0);
        bcast(valid, root);
        if (!valid) return MPI_ERR_OTHER;
        recvbuf.resize(count/nprocs);
        return scatter(sendbuf.data(), count/nprocs, recvbuf.data(), root);
    }

    template <class T> int Comm::scatterv(const std::vector<T>& sendbuf, const std::vector<int>& sendcounts, std::vector<T>& recvbuf, int root) const
    {
        auto comm = getcomm();
        auto dtype = mpi_type<T>();

        int valid = !!(sendcounts.size() == nprocs);
        MPI_Bcast(&valid, 1, MPI_INT, root, comm);

        if (!valid) return MPI_ERR_OTHER;

        int recvcount;
        std::vector<int> displs;

        if (myrank == root)
        {
            displs.resize(nprocs);
            std::exclusive_scan(sendcounts.begin(), sendcounts.end(), displs.begin(), static_cast<int>(0));
        }

        int err = MPI_Scatter(sendcounts.data(), 1, MPI_INT, &recvcount, 1, MPI_INT, root, comm);
        if (err != MPI_SUCCESS) return err;

        recvbuf.resize(recvcount);

        return MPI_Scatterv(sendbuf.data(), sendcounts.data(), displs.data(), dtype, recvbuf.data(), recvcount, dtype, root, comm);
    }

    template <class T> int Comm::scatterv(const std::vector<std::vector<T>>& sendbufs, std::vector<T>& recvbuf, int root) const
    {
        auto comm = getcomm();
        auto dtype = mpi_type<T>();

        int valid = !!(sendbufs.size() == nprocs);
        MPI_Bcast(&valid, 1, MPI_INT, root, comm);

        if (!valid) return MPI_ERR_OTHER;

        int recvcount;
        std::vector<int> sendcounts, displs;
        std::vector<T> sendbuf;

        if (myrank == root)
        {
            sendcounts.resize(nprocs), displs.resize(nprocs);

            for (int i = 0; i < nprocs; ++i)
            {
                sendcounts[i] = sendbufs[i].size();
                displs[i] = (i > 0)? displs[i-1] + sendcounts[i-1] : 0;
            }

            sendbuf.resize(sendcounts.back() + displs.back());

            for (int i = 0; i < nprocs; ++i)
            {
                std::copy(sendbufs[i].begin(), sendbufs[i].end(), sendbuf.begin() + displs[i]);
            }
        }

        int err = MPI_Scatter(sendcounts.data(), 1, MPI_INT, &recvcount, 1, MPI_INT, root, comm);
        if (err != MPI_SUCCESS) return err;

        recvbuf.resize(recvcount);

        return MPI_Scatterv(sendbuf.data(), sendcounts.data(), displs.data(), dtype, recvbuf.data(), recvcount, dtype, root, comm);
    }

    template <class T> int Comm::alltoall(const T* sendbuf, int count, T* recvbuf) const
    {
        return MPI_Alltoall(sendbuf, count, mpi_type<T>(), recvbuf, count, mpi_type<T>(), getcomm());
    }

    template <class T> int Comm::alltoall(const std::vector<T>& sendbuf, std::vector<T>& recvbuf) const
    {
        int valid = !!(sendbuf.size() == nprocs);
        if (!is_same_val(valid)) return MPI_ERR_OTHER;
        if (!valid) return MPI_ERR_OTHER;
        recvbuf.resize(nprocs);
        return alltoall(sendbuf.data(), 1, recvbuf.data());
    }

    template <class T> int Comm::alltoallv(const std::vector<T>& sendbuf, const std::vector<int>& sendcounts, std::vector<T>& recvbuf) const
    {
        auto comm = getcomm();
        auto dtype = mpi_type<T>();

        int valid = !!(sendcounts.size() == nprocs);
        if (!is_same_val(valid)) return MPI_ERR_OTHER;
        if (!valid) return MPI_ERR_OTHER;

        std::vector<int> recvcounts(nprocs), sdispls(nprocs), rdispls(nprocs);

        int err = MPI_Alltoall(sendcounts.data(), 1, MPI_INT, recvcounts.data(), 1, MPI_INT, comm);
        if (err != MPI_SUCCESS) return err;

        std::exclusive_scan(sendcounts.begin(), sendcounts.end(), sdispls.begin(), static_cast<int>(0));
        std::exclusive_scan(recvcounts.begin(), recvcounts.end(), rdispls.begin(), static_cast<int>(0));

        recvbuf.resize(recvcounts.back() + rdispls.back());

        return MPI_Alltoallv(sendbuf.data(), sendcounts.data(), sdispls.data(), dtype,
                             recvbuf.data(), recvcounts.data(), rdispls.data(), dtype, comm);
    }

    template <class T> int Comm::alltoallv(const std::vector<std::vector<T>>& sendbufs, std::vector<T>& recvbuf) const
    {
        auto comm = getcomm();
        auto dtype = mpi_type<T>();

        int valid = !!(sendbufs.size() == nprocs);
        if (!is_same_val(valid)) return MPI_ERR_OTHER;
        if (!valid) return MPI_ERR_OTHER;

        std::vector<T> sendbuf;
        std::vector<int> sendcounts(nprocs), recvcounts(nprocs), sdispls(nprocs), rdispls(nprocs);

        for (int i = 0; i < nprocs; ++i)
        {
            sendcounts[i] = sendbufs[i].size();
            sdispls[i] = (i > 0)? sdispls[i-1] + sendcounts[i-1] : 0;
        }

        int err = MPI_Alltoall(sendcounts.data(), 1, MPI_INT, recvcounts.data(), 1, MPI_INT, comm);
        if (err != MPI_SUCCESS) return err;

        sendbuf.resize(sendcounts.back() + sdispls.back());

        for (int i = 0; i < nprocs; ++i)
        {
            std::copy(sendbufs[i].begin(), sendbufs[i].end(), sendbuf.begin() + sdispls[i]);
        }

        std::exclusive_scan(recvcounts.begin(), recvcounts.end(), rdispls.begin(), static_cast<int>(0));
        recvbuf.resize(recvcounts.back() + rdispls.back());

        return MPI_Alltoallv(sendbuf.data(), sendcounts.data(), sdispls.data(), dtype,
                             recvbuf.data(), recvcounts.data(), rdispls.data(), dtype, comm);
    }

    template <class T> int Comm::balanced_buffer(const std::vector<T>& sendbuf, std::vector<T>& recvbuf) const
    {
        std::vector<std::vector<T>> outgoing(nprocs);

        size_t totsize = sendbuf.size();
        allreduce(totsize, MPI_SUM);

        std::vector<int> counts(nprocs), displs(nprocs);
        get_balanced_counts(counts, totsize);

        std::exclusive_scan(counts.begin(), counts.end(), displs.begin(), static_cast<int>(0));

        size_t mysize = sendbuf.size();
        size_t myoffset;

        exscan(mysize, myoffset, MPI_SUM, (size_t)0);

        for (size_t i = 0; i < mysize; ++i)
        {
            int where = (std::upper_bound(displs.begin(), displs.end(), i+myoffset) - displs.begin()) - 1;
            outgoing[where].push_back(sendbuf[i]);
        }

        return alltoallv(outgoing, recvbuf);
    }

    void Comm::log_strings(const std::string& mystr, std::ostream& os) const
    {
        std::vector<char> mybuf(mystr.begin(), mystr.end());
        std::vector<char> buf;

        gatherv(mybuf, buf, 0);
        std::string str(buf.begin(), buf.end());

        log_string(str, os, 0);
    }

    void Comm::log_string(const std::string& str, std::ostream& os, int root) const
    {
        if (myrank == root)
            os << str;
    }

    template <real_type Real, index_type Index>
    void ArgmaxPair<Real, Index>::mpi_argmax(void *_in, void *_inout, int *len, MPI_Datatype *dtype)
    {
        ArgmaxPair *in = (ArgmaxPair*)_in;
        ArgmaxPair *inout = (ArgmaxPair*)_inout;

        for (int i = 0; i < *len; ++i)
            if (inout[i].value < in[i].value)
            {
                inout[i].value = in[i].value;
                inout[i].index = in[i].index;
            }
    }

    template <real_type Real, index_type Index>
    void ArgmaxPair<Real, Index>::create_mpi_handlers(MPI_Datatype& MPI_ARGMAX_PAIR, MPI_Op& MPI_ARGMAX)
    {
        int blklens[2] = {1,1};
        MPI_Aint disps[2] = {offsetof(ArgmaxPair, index), offsetof(ArgmaxPair, value)};
        MPI_Datatype types[2] = {mpi_type<Index>(), mpi_type<Real>()};
        MPI_Type_create_struct(2, blklens, disps, types, &MPI_ARGMAX_PAIR);
        MPI_Type_commit(&MPI_ARGMAX_PAIR);
        MPI_Op_create(&mpi_argmax, 1, &MPI_ARGMAX);
    }
}
