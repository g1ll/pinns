from datetime import datetime

import numpy as np
import tensorflow as tf
import matplotlib.pyplot as plt

# 1. Configurações do Problema
y0 = 1.0  # Condição inicial: y(0) = 1
k = 0.8   # Constante física da ODE

# Intervalo de tempo para o domínio [0, 2]
t_min, t_max = 0.0, 2.0

# 2. Definição da Arquitetura da Rede Neural
def build_pinn():
    model = tf.keras.Sequential([
        tf.keras.layers.Input(shape=(1,)),
        tf.keras.layers.Dense(32, activation='tanh'),
        tf.keras.layers.Dense(32, activation='tanh'),
        tf.keras.layers.Dense(1, activation=None)
    ])
    return model

pinn = build_pinn()
adam_lrate = 0.0001
optimizer = tf.keras.optimizers.Adam(learning_rate=adam_lrate)

# 3. PONTOS DE COLOCAÇÃO (Collocation Points)
# Pontos internos no domínio [0, 2] onde avaliamos o resíduo da ODE
N_f = 1000
t_f = np.linspace(t_min, t_max, N_f).reshape(-1, 1)
t_f_tf = tf.convert_to_tensor(t_f, dtype=tf.float32)

# Ponto de condição inicial (t = 0)
t_0_tf = tf.constant([[0.0]], dtype=tf.float32)
y_0_tf = tf.constant([[y0]], dtype=tf.float32)

# 4. Função de Perda Física (Physics Loss)
def compute_loss(model, t_f, t_0, y_0):
    # A) Perda da Condição Inicial (Data Loss)
    y_pred_0 = model(t_0)
    loss_init = tf.reduce_mean(tf.square(y_pred_0 - y_0))

    # B) Perda do Resíduo da ODE (Physics Loss)
    with tf.GradientTape() as tape:
        tape.watch(t_f)
        y_pred_f = model(t_f)
    
    # Cálculo automático de dy/dt via TensorFlow Autodiff
    dy_dt = tape.gradient(y_pred_f, t_f)
    
    # Resíduo da física: R(t) = dy/dt + k*y(t) deve ser igual a 0
    residual = dy_dt + k * y_pred_f
    loss_physics = tf.reduce_mean(tf.square(residual))

    # Loss Total
    return loss_init + loss_physics, loss_init, loss_physics

# 5. Passo de treinamento compilado pelo TensorFlow
epochs = 5000
loss_history = []

@tf.function
def train_step():
    with tf.GradientTape() as tape:
        total_loss, loss_init, loss_phys = compute_loss(
            pinn, t_f_tf, t_0_tf, y_0_tf
        )

    grads = tape.gradient(total_loss, pinn.trainable_variables)
    optimizer.apply_gradients(zip(grads, pinn.trainable_variables))
    return total_loss, loss_init, loss_phys


# 6. Loop de Treinamento
for epoch in range(epochs):
    total_loss, loss_init, loss_phys = train_step()

    if epoch % 500 == 0:
        total_loss_value = float(total_loss)
        loss_init_value = float(loss_init)
        loss_phys_value = float(loss_phys)
        print(
            f"Epoch {epoch:4d} | Total Loss: {total_loss_value:.6f} | "
            f"IC Loss: {loss_init_value:.6f} | "
            f"Physics Loss: {loss_phys_value:.6f}"
        )

    if epoch % 10 == 0:
        loss_history.append(float(total_loss))

# 7. Avaliação e Plotagens
t_test = np.linspace(t_min, t_max, 200).reshape(-1, 1)
t_test_tf = tf.convert_to_tensor(t_test, dtype=tf.float32)

# Predição da PINN
y_pinn = pinn(t_test_tf).numpy()

# Solução Analítica Exata
y_exact = y0 * np.exp(-k * t_test)

# Erro Absoluto
abs_error = np.abs(y_exact - y_pinn)

# --- Gráficos de Análise ---
fig, axs = plt.subplots(1, 3, figsize=(16, 4.5))

# Gráfico 1: Comparação da Solução
axs[0].plot(t_test, y_exact, 'r-', linewidth=2, label='Exata: $y_0 e^{-kt}$')
axs[0].plot(t_test, y_pinn, 'k--', linewidth=2, label='PINN')
axs[0].scatter([0], [y0], color='blue', zorder=5, label='Condição Inicial')
axs[0].set_title('Solução da ODE')
axs[0].set_xlabel('Tempo (t)')
axs[0].set_ylabel('y(t)')
axs[0].legend()
axs[0].grid(True, alpha=0.3)

# Gráfico 2: Erro Absoluto
axs[1].plot(t_test, abs_error, 'g-', linewidth=1.5)
axs[1].set_title('Erro Absoluto $|y_{exata} - y_{PINN}|$')
axs[1].set_xlabel('Tempo (t)')
axs[1].set_ylabel('Erro')
axs[1].grid(True, alpha=0.3)

# Gráfico 3: Convergência da Loss
axs[2].semilogy(loss_history, color='purple')
axs[2].set_title(f'Convergência da Loss | Adam LR = {adam_lrate:.4f}')
axs[2].set_xlabel('Época')
axs[2].set_ylabel('Loss (Escala Log)')
axs[2].grid(True, alpha=0.3)

plt.tight_layout()
# plt.show()

timestamp = datetime.now().strftime('%y%m%d_%H%M%S')
output_name = f'resultado_ode_heat_decay_{timestamp}.png'
plt.savefig(output_name, dpi=150, bbox_inches='tight')
print(f"Gráfico de análise salvo com sucesso como '{output_name}'!")