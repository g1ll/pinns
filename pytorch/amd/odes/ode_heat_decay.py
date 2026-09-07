from datetime import datetime

import matplotlib.pyplot as plt
import numpy as np
import torch
from torch import nn


# 1. Configuracoes do problema
y0 = 1.0
k = 0.8
t_min, t_max = 0.0, 2.0

if not torch.cuda.is_available():
    raise RuntimeError(
        "GPU AMD/ROCm nao detectada. Este script nao possui fallback para CPU."
    )
if torch.version.hip is None:
    raise RuntimeError(
        "O PyTorch carregado nao e uma build ROCm/HIP. "
        "Instale o wheel PyTorch para ROCm neste ambiente."
    )

device = torch.device("cuda:0")
print(f"Dispositivo ROCm: {device}")
print(f"GPU: {torch.cuda.get_device_name(device)}")
print(f"ROCm/HIP: {torch.version.hip}")


# 2. Definicao da arquitetura da rede neural
class PINN(nn.Module):
    def __init__(self):
        super().__init__()
        self.network = nn.Sequential(
            nn.Linear(1, 32),
            nn.Tanh(),
            nn.Linear(32, 32),
            nn.Tanh(),
            nn.Linear(32, 1),
        )

    def forward(self, t):
        return self.network(t)


pinn = PINN().to(device)
adam_lrate = 0.0001
optimizer = torch.optim.Adam(pinn.parameters(), lr=adam_lrate)

# 3. Pontos de colocacao no dominio [0, 2]
N_f = 100
t_f = torch.linspace(t_min, t_max, N_f, device=device).reshape(-1, 1)
t_0 = torch.tensor([[0.0]], dtype=torch.float32, device=device)
y_0 = torch.tensor([[y0]], dtype=torch.float32, device=device)


# 4. Funcao de perda fisica
def compute_loss(model, t_f, t_0, y_0):
    # A) Perda da condicao inicial
    y_pred_0 = model(t_0)
    loss_init = torch.mean((y_pred_0 - y_0) ** 2)

    # B) Perda do residuo da ODE
    t_f_for_grad = t_f.detach().clone().requires_grad_(True)
    y_pred_f = model(t_f_for_grad)
    dy_dt = torch.autograd.grad(
        y_pred_f,
        t_f_for_grad,
        grad_outputs=torch.ones_like(y_pred_f),
        create_graph=True,
    )[0]

    # R(t) = dy/dt + k*y(t) deve ser igual a zero
    residual = dy_dt + k * y_pred_f
    loss_physics = torch.mean(residual ** 2)

    return loss_init + loss_physics, loss_init, loss_physics


# 5. Loop de treinamento
epochs = 5000
loss_history = []

for epoch in range(epochs):
    optimizer.zero_grad()
    total_loss, loss_init, loss_phys = compute_loss(pinn, t_f, t_0, y_0)
    total_loss.backward()
    optimizer.step()

    if epoch % 500 == 0:
        print(
            f"Epoch {epoch:4d} | Total Loss {total_loss.item():.6f} | "
            f"IC Loss {loss_init.item():.6f} | "
            f"Physics Loss {loss_phys.item():.6f}"
        )

    if epoch % 10 == 0:
        loss_history.append(total_loss.item())

# 6. Avaliacao e plotagens
t_test = torch.linspace(t_min, t_max, 200, device=device).reshape(-1, 1)
with torch.no_grad():
    y_pinn = pinn(t_test).cpu().numpy()

t_test_np = t_test.cpu().numpy()
y_exact = y0 * np.exp(-k * t_test_np)
abs_error = np.abs(y_exact - y_pinn)

fig, axs = plt.subplots(1, 3, figsize=(16, 4.5))

axs[0].plot(t_test_np, y_exact, "r-", linewidth=2, label="Exata: $y_0 e^{-kt}$")
axs[0].plot(t_test_np, y_pinn, "k--", linewidth=2, label="PINN")
axs[0].scatter([0], [y0], color="blue", zorder=5, label="Condicao Inicial")
axs[0].set_title("Solucao da ODE")
axs[0].set_xlabel("Tempo (t)")
axs[0].set_ylabel("y(t)")
axs[0].legend()
axs[0].grid(True, alpha=0.3)

axs[1].plot(t_test_np, abs_error, "g-", linewidth=1.5)
axs[1].set_title("Erro Absoluto $|y_{exata} - y_{PINN}|$")
axs[1].set_xlabel("Tempo (t)")
axs[1].set_ylabel("Erro")
axs[1].grid(True, alpha=0.3)

axs[2].semilogy(loss_history, color="purple")
axs[2].set_title(f"Convergencia da Loss | Adam LR = {adam_lrate:.4f}")
axs[2].set_xlabel("Epoca")
axs[2].set_ylabel("Loss (Escala Log)")
axs[2].grid(True, alpha=0.3)

plt.tight_layout()
timestamp = datetime.now().strftime("%y%m%d_%H%M%S")
output_name = f"./odes/results/resultado_ode_heat_decay_{timestamp}.png"
plt.savefig(output_name, dpi=150, bbox_inches="tight")
print(f"Grafico de analise salvo com sucesso como '{output_name}'!")