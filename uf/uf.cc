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
		int find_state_flag = 0;
		int find_parent_flag = 0;
		int find_rsp_state_flag = 0;
		int find_rsp_parent_flag = 0;
		int find_state_tag = 0;
		int find_rsp_state_tag = 512;
		int find_parent_tag = 1024;
		int find_rsp_parent_tag = 2048;
		MPI_Status state_status;
		MPI_Status parent_status;
		MPI_Message state_handle;
		MPI_Message parent_handle;
		MPI_Status state_rsp_status;
		MPI_Status parent_rsp_status;
		MPI_Message state_rsp_handle;
		MPI_Message parent_rsp_handle;

		MPI_Comm_size(MPI_COMM_WORLD, &size);
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);

		for (int i = 0; i < size; i++) {
			MPI_Improbe(i, find_state_tag + i, MPI_COMM_WORLD, &find_state_flag,
					&state_handle, &state_status);
			MPI_Improbe(i, find_parent_tag + i, MPI_COMM_WORLD,
					&find_parent_flag, &parent_handle, &parent_status);

			if (find_state_flag) {
				int state_message;
				int state;

				MPI_Mrecv(&state_message, 1, MPI_INT, &state_handle,
						&state_status);
				state = state_message;
				this->validate(state);

				while (state != parent[state]) {
					if (parent[state] % size != rank) {
						int parent_message;
						MPI_Status recv_status;

						MPI_Send(&parent[state], 1, MPI_INT, state % size,
								find_parent_tag + rank, MPI_COMM_WORLD);
						MPI_Recv(&parent_message, 1, MPI_INT, state % size,
								find_rsp_parent_tag + state % size,
								MPI_COMM_WORLD, &recv_status);
						parent[state] = parent_message;
						state = parent[state];
						continue;
					}

					parent[state] = parent[parent[state]];
					state = parent[state];
				}

				MPI_Send(&state, 1, MPI_INT, i, find_rsp_state_tag + i,
				MPI_COMM_WORLD);
			}

			if (find_parent_flag) {
				int parent_message;
				int parent;

				MPI_Mrecv(&parent_message, 1, MPI_INT, &parent_handle,
						&parent_status);
				parent = parent[parent_message];
				MPI_Send(&parent, 1, MPI_INT, i, find_rsp_parent_tag + i,
				MPI_COMM_WORLD);
			}
		}

		if (p % size != rank) {
			int state_message;
			MPI_Request request;

			MPI_Isend(&p, 1, MPI_INT, p % size, find_state_tag + rank,
			MPI_COMM_WORLD, &request);

			while (!find_rsp_state_flag) {
				MPI_Improbe(p % size, find_rsp_state_tag + p % size,
				MPI_COMM_WORLD, &find_rsp_state_flag, &state_rsp_handle,
						&state_rsp_status);

				for (int i = 0; i < size; i++) {
					MPI_Improbe(i, find_rsp_parent_tag + i, MPI_COMM_WORLD,
							&find_rsp_parent_flag, &parent_rsp_handle,
							&parent_rsp_status);

					if (find_rsp_parent_flag) {
						int parent_message;
						int parent;

						MPI_Mrecv(&parent_message, 1, MPI_INT, &parent_handle,
								&parent_status);
						parent = parent[parent_message];
						MPI_Send(&parent, 1, MPI_INT, i,
								find_rsp_parent_tag + i, MPI_COMM_WORLD); // i
					}
				}

				if (find_rsp_state_flag) {
					MPI_Mrecv(&state_message, 1, MPI_INT, &state_rsp_handle,
							&state_rsp_status);
				}
			}

			p = state_message;
			return p;
		}

		this->validate(p);

		while (p != parent[p]) {
			if (parent[p] % size != rank) {
				int parent_message;
				MPI_Status recv_status;

				MPI_Send(&parent[p], 1, MPI_INT, parent[p] % size,
						find_parent_tag + rank, MPI_COMM_WORLD);
				MPI_Recv(&parent_message, 1, MPI_INT, parent[p] % size,
						find_rsp_parent_tag + parent[p] % size, MPI_COMM_WORLD,
						&recv_status);
				parent[p] = parent_message;
				p = parent[p];
				continue;
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

