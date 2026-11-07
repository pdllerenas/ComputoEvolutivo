import numpy as np
import matplotlib.pyplot as plt

# 1. Definir la malla de puntos para el espacio de búsqueda (x >= 0, y >= 0)
x = np.linspace(0, 12, 500)
y = np.linspace(0, 12, 500)
X, Y = np.meshgrid(x, y)

# 2. Definir las restricciones en formato h(x,y) <= 0
H1 = (X - 5)**2 + Y**2 - 26
H2 = 4*X + Y - 20

# 3. Identificar la región factible
factible = (H1 <= 0) & (H2 <= 0) & (X >= 0) & (Y >= 0)

plt.figure(figsize=(10, 8))

# 4. Sombrear la región factible
plt.imshow(factible, extent=(0, 12, 0, 12), origin='lower', alpha=0.3, cmap='Greens')

# 5. Dibujar contornos de las fronteras de las restricciones
c1 = plt.contour(X, Y, H1, levels=[0], colors='blue', linewidths=2)
c2 = plt.contour(X, Y, H2, levels=[0], colors='red', linewidths=2)
plt.clabel(c1, fmt='h1(x,y)=0', inline=True, fontsize=11)
plt.clabel(c2, fmt='h2(x,y)=0', inline=True, fontsize=11)

# 6. Datos de los puntos calculados a mano
puntos_orig = [(10, 5), (5, 5), (2, 8), (8, 8)]
puntos_nuevos = [(2.467, 10.133), (3.824, 4.706), (2.966, 5.425), (3.293, 6.828)]
labels = ['P1', 'P2', 'P3', 'P4']

# 7. Graficar puntos y sus vectores de reparación (Takahama)
for (po, pn, label) in zip(puntos_orig, puntos_nuevos, labels):
    plt.plot(po[0], po[1], 'ko', markersize=6) # Punto original infactible
    plt.plot(pn[0], pn[1], 'b*', markersize=8) # Punto reparado
    plt.text(po[0]+0.2, po[1]+0.2, label, fontsize=12, fontweight='bold')
    plt.annotate('', xy=pn, xytext=po, 
                 arrowprops=dict(arrowstyle="->", color='black', lw=1.5))

# 8. Marcar el mínimo global conocido
plt.plot(3, 2, 'rP', markersize=12, label='Mínimo Óptimo (3,2)')

# 9. Configuraciones estéticas de la gráfica
plt.xlim(0, 12)
plt.ylim(0, 12)
plt.xlabel('Eje X', fontsize=12)
plt.ylabel('Eje Y', fontsize=12)
plt.title('Takahama\nFeasible region (green), infeaseible (white)', fontsize=14)
plt.grid(True, linestyle='--', alpha=0.6)
plt.legend(loc='upper right')
plt.tight_layout()

# Mostrar gráfica
plt.savefig("takahama.png")
plt.show()
