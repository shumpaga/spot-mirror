#include <iostream>
#include "mpi.h"

class union_find {
public:

	// n, number of entries
	union_find(int n) {
		if (n < 0)
			this->exit_code = 1;

		this->count = n;
		this->size = n;
		this->parent = new int[n];
		this->rank = new int[n];

		for (int i = 0; i < n; i++) {
			this->parent[i] = i;
			this->rank[i] = 0;
		}
	}

	~union_find(void) {
		delete[] this->parent;
		delete[] this->rank;
	}

	// return, the component identifier for the component containing the entry p
	int find(int p) {
		int size;
		int rank;
		int flag = 0;
		int find_tag = 0;
		int find_rsp_tag = 1;
		MPI_Status send_status;
		MPI_Message handle;

		MPI_Comm_size(MPI_COMM_WORLD, &size);
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);

		for (int i = 0; i < size; i++) {
			MPI_Improbe(i, find_tag + i, MPI_COMM_WORLD, &flag, &handle,
					&send_status);

			if (flag) {
				int message;

				MPI_Mrecv(&message, 1, MPI_INT, &handle, &send_status);

				this->validate(message);

				while (message != parent[message]) {
					parent[message] = parent[parent[message]];
					message = parent[message];
				}

				MPI_Send(&message, 1, MPI_INT, i, find_rsp_tag, MPI_COMM_WORLD);
			}
		}

		if (p % size != rank) {
			int message;

			MPI_Status recv_status;
			MPI_Send(&p, 1, MPI_INT, p % size, find_tag, MPI_COMM_WORLD);
			MPI_Recv(&message, 1, MPI_INT, p % size, find_rsp_tag,
			MPI_COMM_WORLD, &recv_status);
			p = message;
		}

		this->validate(p);

		while (p != parent[p]) {
			if (parent[p] % size != rank) {
				int message;

				MPI_Status recv_status;
				MPI_Send(&p, 1, MPI_INT, p % size, find_tag, MPI_COMM_WORLD);
				MPI_Recv(&message, 1, MPI_INT, p % size, find_rsp_tag,
				MPI_COMM_WORLD, &recv_status);
				parent[p] = message;
			}

			parent[p] = parent[parent[p]];
			p = parent[p];
		}

		return p;
	}

	// return, number of components between 1 and size
	int count_components(void) {
		return this->count;
	}

	// return, true if the the two entries are in the same component
	bool connected(int p, int q) {
		return find(p) == find(q);
	}

	//  Merges the component containing p with the the component containing q
	void unite(int p, int q) {
		int root_p = find(p);
		int root_q = find(q);

		if (root_p == root_q)
			return;

		if (this->rank[root_p] < this->rank[root_q])
			this->parent[root_p] = root_q;

		else {
			if (this->rank[root_q] < this->rank[root_p])
				this->parent[root_q] = root_p;

			else {
				this->parent[root_q] = root_p;
				this->rank[root_p]++;
			}
		}

		this->count--;
		return;
	}

private:
	int size;
	int exit_code;
	int count;
	int* parent;
	int * rank;

	// validate that p is a valid index
	void validate(int p) {
		int n = this->size;

		if (p < 0 || p >= n)
			std::cerr << "index " + p + " is not between 0 and " + (n - 1)
					<< std::endl;
	}
};

int main(int argc, char** argv) {
	union_find uf;
	int n;
	int p;
	int q;

	MPI_Init(&argc, &argv);

	std::cin >> n;
	uf = new union_find(n);

	while (std::cin >> p >> q) {
		if (uf.connected(p, q))
			continue;

		uf.unite(p, q);
		std::cout << p << " " << q << std::endl;
	}

	std::cout << uf.count_components() << "components" << std::endl;
	delete uf;
	return 0;
}

