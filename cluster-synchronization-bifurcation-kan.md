# Cluster Synchronization, Bifurcation Theory, and KAN
# Đồng bộ hóa cụm, lý thuyết điểm rẽ nhánh và KAN

## Scope and important qualifications / Phạm vi và lưu ý quan trọng

**English.** This note gives a concrete four-node example and an executable learning prototype. A Fiedler eigenvalue measures graph connectivity; by itself it does **not** determine the synchronization threshold of a Kuramoto system. The threshold also depends on the oscillator frequencies, edge weights, and coupling law. Likewise, finite time-series data do not uniquely identify an arbitrary interaction law without assumptions and sufficiently informative observations. A KAN can approximate a learned law, but symbolic recovery is not guaranteed by the Kolmogorov–Arnold representation theorem.

**Tiếng Việt.** Tài liệu này trình bày ví dụ cụ thể với bốn nút và một nguyên mẫu học có thể chạy. Trị riêng Fiedler đo độ kết nối của đồ thị; riêng nó **không** xác định ngưỡng đồng bộ của hệ Kuramoto. Ngưỡng còn phụ thuộc vào tần số dao động, trọng số cạnh và quy luật ghép nối. Tương tự, dữ liệu chuỗi thời gian hữu hạn không xác định duy nhất một hàm tương tác bất kỳ nếu thiếu giả định và dữ liệu quan sát đủ thông tin. KAN có thể xấp xỉ quy luật đã học, nhưng định lý biểu diễn Kolmogorov–Arnold không đảm bảo tự khôi phục được công thức ký hiệu chính xác.

## 1. Model / Mô hình

For an undirected weighted graph, consider the sine-coupled Kuramoto model

$$
\dot\theta_i = \omega_i + K\sum_{j=1}^N A_{ij}\sin(\theta_j-\theta_i),\qquad i=1,\ldots,N.
$$

Here, $\theta_i$ is the phase, $\omega_i$ its natural frequency, $A$ the symmetric nonnegative adjacency matrix, and $K\geq0$ the global coupling strength. A phase-locked state has constant phase differences and a common rotation frequency; it need not have equal phases. The graph Laplacian is $L=D-A$, where $D_{ii}=\sum_j A_{ij}$.

Xét mô hình Kuramoto ghép nối bằng hàm sin trên đồ thị vô hướng có trọng số:

$$
\dot\theta_i = \omega_i + K\sum_{j=1}^N A_{ij}\sin(\theta_j-\theta_i),\qquad i=1,\ldots,N.
$$

Trong đó $\theta_i$ là pha, $\omega_i$ là tần số tự nhiên, $A$ là ma trận kề đối xứng không âm, và $K\geq0$ là cường độ ghép nối toàn cục. Trạng thái khóa pha có hiệu pha không đổi và cùng tần số quay; các pha không nhất thiết bằng nhau. Ma trận Laplace là $L=D-A$, với $D_{ii}=\sum_j A_{ij}$.

## 2. Four-node example and threshold / Ví dụ bốn nút và ngưỡng

Take two internally coupled pairs, with unit weights on edges $(1,2)$ and $(3,4)$ and a bridge of weight $\epsilon=0.1$ on $(2,3)$:

$$
A=\begin{pmatrix}0&1&0&0\\1&0&0.1&0\\0&0.1&0&1\\0&0&1&0\end{pmatrix},\qquad
L=\begin{pmatrix}1&-1&0&0\\-1&1.1&-0.1&0\\0&-0.1&1.1&-1\\0&0&-1&1\end{pmatrix}.
$$

The Laplacian eigenvalues are approximately $0, 0.095012, 2, 2.104988$, so the algebraic connectivity is $\lambda_2\approx0.095012$. This indicates a weakly connected bottleneck, not a universal critical-coupling formula.

For the specific frequencies $\omega=(0.19,0.19,-0.19,-0.19)$, the sum is zero. In a locked state the required steady flow across each successive edge is determined by the cumulative frequency on one side of that edge: magnitudes $0.19$, $0.38$, and $0.19$. An edge of weight $w$ can carry at most $Kw$, since $|\sin(\Delta\theta)|\leq1$. Thus the exact locking condition for this tree and this frequency vector is

$$
K > \max\left(\frac{0.19}{1},\frac{0.38}{0.1},\frac{0.19}{1}\right)=3.8.
$$

At $K=3.8$ the bridge is at its limiting flow ($|\sin(\Delta\theta)|=1$); a stable locked solution exists for $K>3.8$. This result follows from flow balance on this **tree** and the sine-coupling capacity. It is not obtained by substituting $\lambda_2$ into $\Delta\omega/(2\lambda_2)$. For general graphs, locking thresholds require analysis of the network flow constraints and stability; the Fiedler value alone is insufficient.

Với tần số cụ thể $\omega=(0.19,0.19,-0.19,-0.19)$, tổng tần số bằng không. Trong trạng thái khóa, dòng ghép nối ổn định qua từng cạnh liên tiếp được xác định bởi tổng tần số ở một phía cạnh đó: độ lớn lần lượt là $0.19$, $0.38$ và $0.19$. Cạnh trọng số $w$ truyền được nhiều nhất $Kw$, vì $|\sin(\Delta\theta)|\leq1$. Do đó, điều kiện khóa chính xác cho cây này và véc-tơ tần số này là

$$
K > \max\left(\frac{0.19}{1},\frac{0.38}{0.1},\frac{0.19}{1}\right)=3.8.
$$

Tại $K=3.8$, cầu nối đã đạt giới hạn dòng ($|\sin(\Delta\theta)|=1$); trạng thái khóa ổn định tồn tại khi $K>3.8$. Kết quả này dựa trên cân bằng dòng trên **đồ thị cây** và giới hạn của ghép nối sin. Không thể suy ra nó bằng cách thế $\lambda_2$ vào công thức $\Delta\omega/(2\lambda_2)$. Với đồ thị tổng quát, cần phân tích các ràng buộc dòng mạng và độ ổn định; riêng trị riêng Fiedler là chưa đủ.

## 3. Graph-KAN and identifiability / Graph-KAN và tính khả định

A graph model can assign a learnable univariate function $\Phi_{ij}$ to each candidate edge and predict $\dot\theta_i=\omega_i+\sum_j\Phi_{ij}(\theta_j-\theta_i)$. A sparsity penalty may encourage unused candidate edges to shrink. For undirected sine coupling, the interaction should respect the antisymmetry $\Phi_{ji}(-x)=-\Phi_{ij}(x)$; an unconstrained model may not preserve this property.

A practical workflow is to estimate derivatives from sufficiently sampled trajectories (or fit an ODE solver directly), train on multiple initial conditions, validate on held-out trajectories, and check sensitivity to noise and regularization. A spline or radial-basis fit is an approximation. Matching it to a symbolic library is a separate model-selection step; report the candidate library, fit error, and uncertainty. Recovering $K A_{ij}\sin(x)$ exactly is only justified when the data and assumptions identify that law.

Mô hình đồ thị có thể gán một hàm đơn biến khả học $\Phi_{ij}$ cho mỗi cạnh ứng viên và dự đoán $\dot\theta_i=\omega_i+\sum_j\Phi_{ij}(\theta_j-\theta_i)$. Phạt thưa có thể khuyến khích các cạnh không dùng co về không. Với ghép nối sin vô hướng, tương tác cần thỏa tính phản đối xứng $\Phi_{ji}(-x)=-\Phi_{ij}(x)$; mô hình không ràng buộc có thể không bảo toàn tính chất này.

Quy trình thực tế là ước lượng đạo hàm từ quỹ đạo được lấy mẫu đủ dày (hoặc khớp trực tiếp bộ giải ODE), huấn luyện trên nhiều điều kiện đầu, kiểm định trên quỹ đạo giữ lại và kiểm tra độ nhạy với nhiễu cũng như hệ số điều chuẩn. Khớp spline hoặc cơ sở xuyên tâm chỉ là xấp xỉ. Đối chiếu với thư viện ký hiệu là bước chọn mô hình riêng; cần báo cáo thư viện ứng viên, sai số khớp và độ bất định. Chỉ có thể kết luận khôi phục chính xác $K A_{ij}\sin(x)$ khi dữ liệu và giả định đủ để xác định quy luật đó.

## 4. Executable reference implementation / Mã tham khảo có thể thực thi

The script below computes the spectrum and the exact tree threshold for this example, then simulates the system above threshold. It deliberately does not claim that a KAN can infer an exact law from one trajectory.

Đoạn mã dưới đây tính phổ và ngưỡng chính xác của cây trong ví dụ, sau đó mô phỏng hệ ở phía trên ngưỡng. Mã không khẳng định KAN có thể suy ra chính xác quy luật chỉ từ một quỹ đạo.

```python
import numpy as np
from scipy.integrate import solve_ivp


def analyze_cluster_bifurcation(epsilon=0.1, delta_omega=0.38):
    """Return adjacency, Laplacian spectrum, Fiedler value, and tree threshold."""
    if epsilon <= 0:
        raise ValueError("epsilon must be positive for a connected bridge")
    if delta_omega < 0:
        raise ValueError("delta_omega must be nonnegative")

    A = np.array([
        [0.0, 1.0, 0.0, 0.0],
        [1.0, 0.0, epsilon, 0.0],
        [0.0, epsilon, 0.0, 1.0],
        [0.0, 0.0, 1.0, 0.0],
    ])
    L = np.diag(A.sum(axis=1)) - A
    eigenvalues = np.linalg.eigvalsh(L)
    lambda_2 = eigenvalues[1]

    # For omega=(d/2, d/2, -d/2, -d/2), the bridge carries total flow d.
    # The intra-cluster edges each carry d/2; all edge capacities must hold.
    K_critical = max(delta_omega / 2.0, delta_omega / epsilon)
    return A, eigenvalues, lambda_2, K_critical


def main():
    A, eigenvalues, lambda_2, K_critical = analyze_cluster_bifurcation()
    print("Eigenvalues:", np.round(eigenvalues, 6))
    print(f"Fiedler value: {lambda_2:.6f}")
    print(f"Tree locking boundary: K = {K_critical:.6f}")

    omega = np.array([0.19, 0.19, -0.19, -0.19])
    K = 4.5  # Above the exact threshold 3.8 for this example.

    def kuramoto_rhs(_t, theta):
        phase_difference = theta[None, :] - theta[:, None]
        return omega + K * np.sum(A * np.sin(phase_difference), axis=1)

    times = np.linspace(0.0, 20.0, 1001)
    solution = solve_ivp(
        kuramoto_rhs,
        (times[0], times[-1]),
        [0.0, 0.2, 1.5, 1.7],
        t_eval=times,
        rtol=1e-9,
        atol=1e-11,
    )
    if not solution.success:
        raise RuntimeError(solution.message)

    print("Simulated phase array shape:", solution.y.T.shape)
    print("Final phase differences:", np.round(
        np.diff(solution.y[:, -1]), 6
    ))


if __name__ == "__main__":
    main()
```

## 5. References and interpretation / Tài liệu tham khảo và diễn giải

The Kuramoto model is a standard phase-oscillator model; the graph Laplacian and its second-smallest eigenvalue are standard tools in spectral graph theory. Neither the Russian mathematical tradition nor the Kolmogorov–Arnold theorem implies the specific threshold or automatic discovery claims corrected above. The numerical and analytical statements in this note are scoped to the stated graph, coupling law, and frequency vector.

Mô hình Kuramoto là mô hình chuẩn của các dao động tử pha; ma trận Laplace đồ thị và trị riêng dương nhỏ nhất là các công cụ chuẩn của lý thuyết phổ đồ thị. Trường phái toán học Nga hay định lý Kolmogorov–Arnold không tự suy ra ngưỡng cụ thể hoặc khả năng khám phá tự động như các tuyên bố đã được hiệu chỉnh ở trên. Các kết quả số và giải tích trong tài liệu này chỉ áp dụng cho đồ thị, quy luật ghép nối và véc-tơ tần số đã nêu.
