import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
import numpy as np
import sys

# ==========================================
# 1. CONFIGURACION Y ESTILOS
# ==========================================
TEXTO_PIE_PAGINA = (
    "Fuentes: Ecuaciones Hyland & Wexler (1983) - ASHRAE Handbook Fundamentals.\n"
    "Desarrollado con Python.\n"
    "Autor: Marcos Sánchez."
)

COLOR_SATURACION = 'black'
COLOR_ENTALPIA = '#FF8C00'  # Naranja (Eje y lineas)
COLOR_HR = '#A52A2A'  # Marron (Humedad Relativa)
COLOR_TBH = 'green'    # Verde (Bulbo Humedo)
COLOR_VOL = 'blue'     # Azul (Volumen Especifico)
COLOR_GRILLA = 'gray'     # Grillado tenue
COLOR_PUNTO = 'red'

PASO_ENTALPIA = 5
PASO_TBH = 5
PASO_VOLUMEN = 0.05

OPCIONES_GUARDADO = {
    "Pequeno (Web)": (10, 8, 100),
    "Mediano (A4)": (11.69, 8.27, 300),
    "Poster":  (24, 18, 300)
}

# ==========================================
# 2. FORMULAS FISICAS (ASHRAE HYLAND & WEXLER)
# ==========================================


def presion_saturacion_ashrae(t_c):
    """
    Calcula Psat usando coeficientes de Hyland & Wexler (1983).
    Estandar ASHRAE Handbook.
    """
    tk = np.array(t_c) + 273.15

    # Coeficientes Agua Liquida (0 a 200 C)
    C8, C9, C10 = -5.8002206e03, 1.3914993, -4.8640239e-02
    C11, C12, C13 = 4.1764768e-05, -1.4452093e-08, 6.5459673

    # Coeficientes Hielo (-100 a 0 C)
    C1, C2, C3 = -5.6745359e03, 6.3925247, -9.6778430e-03
    C4, C5, C6, C7 = 6.2215701e-07, 2.0747825e-09, -9.4840240e-13, 4.1635019

    log_p_ice = (C1/tk + C2 + C3*tk + C4*tk**2 +
                 C5*tk**3 + C6*tk**4 + C7*np.log(tk))
    p_ice = np.exp(log_p_ice)

    log_p_liq = (C8/tk + C9 + C10*tk + C11*tk**2 + C12*tk**3 + C13*np.log(tk))
    p_liq = np.exp(log_p_liq)

    if np.ndim(tk) == 0:
        return p_liq if t_c >= 0 else p_ice
    return np.where(t_c >= 0, p_liq, p_ice)


def presion_vapor_parcial(t_c, hr):
    return hr * presion_saturacion_ashrae(t_c)


def razon_humedad_desde_pv(pv, presion):
    denom = presion - pv
    return 0.621945 * pv / np.maximum(denom, 1.0)


def razon_humedad_desde_hr(t_c, hr, presion):
    es = presion_saturacion_ashrae(t_c)
    pw = hr * es
    return razon_humedad_desde_pv(pw, presion)


def razon_humedad_saturacion(t_c, presion):
    return razon_humedad_desde_hr(t_c, 1.0, presion)


def entalpia_aire(t_c, W):
    return 1.006 * t_c + W * (2501.0 + 1.86 * t_c)


def volumen_especifico(t_c, W, presion):
    return (287.055 * (t_c + 273.15) / presion) * (1 + 1.6078 * W)


def temperatura_rocio(pv):
    # Estimacion precisa inversa usando Magnus para visualizacion rapida
    pv_safe = np.maximum(pv, 1e-6)
    y = np.log(pv_safe / 610.94)
    return 243.04 * y / (17.625 - y)


def temperatura_bulbo_humedo(t_c, hr, presion):
    w_actual = razon_humedad_desde_hr(t_c, hr, presion)
    tw_min, tw_max = -50, t_c
    for _ in range(25):
        tw_mid = (tw_min + tw_max) / 2
        w_s_mid = razon_humedad_saturacion(tw_mid, presion)
        h_s_mid = entalpia_aire(tw_mid, w_s_mid)
        w_calc = (h_s_mid - 1.006 * t_c) / (2501.0 + 1.86 * t_c)
        if w_calc < w_actual:
            tw_min = tw_mid
        else:
            tw_max = tw_mid
    return tw_max


def calor_especifico_aire(W):
    return 1.006 + 1.86 * W


def resolver_t_para_h_y_w(h, w_objetivo):
    denom = 1.006 + 1.86 * w_objetivo
    if denom == 0:
        return 0
    return (h - 2501.0 * w_objetivo) / denom


def obtener_t_saturacion_en_w(w_objetivo, presion):
    t_min, t_max = -50, 150
    for _ in range(50):
        t_med = (t_min + t_max) / 2
        w_val = razon_humedad_saturacion(t_med, presion)
        if w_val < w_objetivo:
            t_min = t_med
        else:
            t_max = t_med
    return t_med


def temperatura_desde_entalpia_saturada(h, presion):
    t_min, t_max = -50, 150
    for _ in range(50):
        t_med = (t_min + t_max) / 2
        w_med = razon_humedad_saturacion(t_med, presion)
        h_med = entalpia_aire(t_med, w_med)
        if h_med < h:
            t_min = t_med
        else:
            t_max = t_med
    return t_min

# ==========================================
# 3. CLASE CALCULADORA
# ==========================================


class CalculadoraPunto:
    def __init__(self, parent, presion_pa, callback_graficar):
        self.ventana = tk.Toplevel(parent)
        self.ventana.title("Calculadora ASHRAE")
        self.ventana.geometry("450x650")
        self.presion = presion_pa
        self.callback = callback_graficar

        frame_modo = ttk.LabelFrame(
            self.ventana, text="Datos de Entrada", padding=10)
        frame_modo.pack(fill='x', padx=10, pady=5)

        self.modo_var = tk.StringVar(value="TBS_HR")
        modos = [
            ("T. Bulbo Seco + Humedad Relativa", "TBS_HR"),
            ("T. Bulbo Seco + T. Bulbo Humedo", "TBS_TBH"),
            ("T. Bulbo Seco + T. Rocio", "TBS_TR"),
            ("T. Bulbo Seco + Humedad Absoluta", "TBS_W")
        ]
        for texto, valor in modos:
            ttk.Radiobutton(frame_modo, text=texto, variable=self.modo_var,
                            value=valor, command=self.actualizar_inputs).pack(anchor='w')

        frame_inputs = ttk.Frame(self.ventana, padding=10)
        frame_inputs.pack(fill='x', padx=10)
        self.lbl_1 = ttk.Label(frame_inputs, text="T. Bulbo Seco (°C):")
        self.lbl_1.grid(row=0, column=0, sticky='w')
        self.ent_1 = ttk.Entry(frame_inputs)
        self.ent_1.grid(row=0, column=1, sticky='ew', padx=5)
        self.lbl_2 = ttk.Label(frame_inputs, text="Humedad Relativa (%):")
        self.lbl_2.grid(row=1, column=0, sticky='w')
        self.ent_2 = ttk.Entry(frame_inputs)
        self.ent_2.grid(row=1, column=1, sticky='ew', padx=5)

        btn_calc = ttk.Button(
            frame_inputs, text="Calcular", command=self.calcular)
        btn_calc.grid(row=2, column=0, columnspan=2, pady=10, sticky='ew')

        self.tree = ttk.Treeview(self.ventana, columns=(
            "Valor", "Unidad"), show="tree headings")
        self.tree.heading("#0", text="Propiedad")
        self.tree.heading("Valor", text="Valor")
        self.tree.heading("Unidad", text="Unidad")
        self.tree.column("#0", width=180)
        self.tree.pack(fill='both', expand=True, padx=10, pady=5)
        self.actualizar_inputs()

    def actualizar_inputs(self):
        modo = self.modo_var.get()
        self.lbl_1.config(text="T. Bulbo Seco (°C):")
        self.ent_1.delete(0, tk.END)
        self.ent_2.delete(0, tk.END)
        if modo == "TBS_HR":
            self.lbl_2.config(text="Humedad Relativa (%):")
        elif modo == "TBS_TBH":
            self.lbl_2.config(text="T. Bulbo Humedo (°C):")
        elif modo == "TBS_TR":
            self.lbl_2.config(text="T. Rocio (°C):")
        elif modo == "TBS_W":
            self.lbl_2.config(text="Humedad Abs. (g/kg):")

    def calcular(self):
        try:
            if not self.ent_1.get() or not self.ent_2.get():
                return
            v1, v2 = float(self.ent_1.get()), float(self.ent_2.get())
            modo = self.modo_var.get()
            tbs, hr, w = v1, 0, 0

            if modo == "TBS_HR":
                hr = v2 / 100.0
                if not (0 <= hr <= 1.0):
                    raise ValueError("HR entre 0 y 100")
                w = razon_humedad_desde_hr(tbs, hr, self.presion)
            elif modo == "TBS_TBH":
                tbh = v2
                if tbh > tbs:
                    raise ValueError("TBH > TBS")
                es_wb = presion_saturacion_ashrae(tbh)
                denom = self.presion - es_wb
                w_sat_wb = 0.621945 * es_wb / np.maximum(denom, 1.0)
                h_sat_wb = entalpia_aire(tbh, w_sat_wb)
                w = (h_sat_wb - 1.006 * tbs) / (2501.0 + 1.86 * tbs)
                pv = w * self.presion / (0.621945 + w)
                es = presion_saturacion_ashrae(tbs)
                hr = pv / es
            elif modo == "TBS_TR":
                pv = presion_saturacion_ashrae(v2)
                w = razon_humedad_desde_pv(pv, self.presion)
                es = presion_saturacion_ashrae(tbs)
                hr = pv / es
            elif modo == "TBS_W":
                w = v2 / 1000.0
                es = presion_saturacion_ashrae(tbs)
                pv = w * self.presion / (0.621945 + w)
                hr = pv / es

            pv = hr * presion_saturacion_ashrae(tbs)
            h = entalpia_aire(tbs, w)
            vol = volumen_especifico(tbs, w, self.presion)
            t_rocio = temperatura_rocio(pv)
            tbh_val = temperatura_bulbo_humedo(
                tbs, hr, self.presion) if modo != "TBS_TBH" else v2
            densidad = 1 / vol if vol > 0 else 0
            p_sat_tbs = presion_saturacion_ashrae(tbs)
            cp = calor_especifico_aire(w)

            self.tree.delete(*self.tree.get_children())
            datos = [
                ("Presion Atmosferica", f"{self.presion:.0f}", "Pa"),
                ("T. Bulbo Seco", f"{tbs:.2f}", "°C"),
                ("Humedad Absoluta", f"{w*1000:.2f}", "g/kg"),
                ("Humedad Relativa", f"{hr*100:.1f}", "%"),
                ("T. Bulbo Humedo", f"{tbh_val:.2f}", "°C"),
                ("T. Rocio", f"{t_rocio:.2f}", "°C"),
                ("Entalpia", f"{h:.2f}", "kJ/kg"),
                ("Presion Vapor", f"{pv:.0f}", "Pa"),
                ("Presion Sat. Vapor", f"{p_sat_tbs:.0f}", "Pa"),
                ("Volumen Especifico", f"{vol:.3f}", "m³/kg"),
                ("Calor Especifico", f"{cp:.3f}", "kJ/(kg.K)"),
                ("Densidad", f"{densidad:.3f}", "kg/m³")
            ]
            for p, v, u in datos:
                self.tree.insert("", tk.END, text=p, values=(v, u))
            self.callback(tbs, w*1000.0)
        except ValueError as e:
            messagebox.showerror("Error", str(e))
        except Exception as e:
            messagebox.showerror("Error", f"Entrada invalida: {e}")

# ==========================================
# 4. APLICACION PRINCIPAL
# ==========================================


class AplicacionPsicrometrica:
    def __init__(self, raiz):
        self.raiz = raiz
        self.raiz.title("Carta Psicrometrica Pro (Hyland & Wexler)")
        self.raiz.protocol("WM_DELETE_WINDOW", self.cerrar_aplicacion)

        self.punto_actual = None
        self.def_altitud = 124
        self.def_t_min = -10
        self.def_t_max = 55
        self.def_w_max = 30

        marco_control = ttk.Frame(raiz, padding="10")
        marco_control.pack(side=tk.LEFT, fill=tk.Y)
        ttk.Label(marco_control, text="Configuracion",
                  font=('Arial', 12, 'bold')).pack(pady=5)
        self.ent_alt = self.crear_input(
            marco_control, "Altitud (msnm):", self.def_altitud)
        self.ent_tmin = self.crear_input(
            marco_control, "T. Minima (C):", self.def_t_min)
        self.ent_tmax = self.crear_input(
            marco_control, "T. Maxima (C):", self.def_t_max)
        self.ent_wmax = self.crear_input(
            marco_control, "W Max (g/kg):", self.def_w_max)

        btn_act = ttk.Button(
            marco_control, text="Actualizar Grafico", command=self.actualizar_grafico)
        btn_act.pack(pady=15, fill=tk.X)
        ttk.Separator(marco_control, orient='horizontal').pack(
            fill='x', pady=10)

        ttk.Label(marco_control, text="Herramientas",
                  font=('Arial', 12, 'bold')).pack(pady=5)
        ttk.Button(marco_control, text="Calculadora de Punto",
                   command=self.abrir_calculadora).pack(pady=5, fill=tk.X)

        ttk.Label(marco_control, text="Exportar:").pack(
            anchor='w', pady=(10, 0))
        self.combo_tamano = ttk.Combobox(marco_control, values=list(
            OPCIONES_GUARDADO.keys()), state="readonly")
        self.combo_tamano.current(0)
        self.combo_tamano.pack(fill=tk.X, pady=5)
        ttk.Button(marco_control, text="Guardar Imagen",
                   command=self.guardar_grafico).pack(pady=5, fill=tk.X)

        marco_plot = ttk.Frame(raiz)
        marco_plot.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True)

        self.figura, self.ejes = plt.subplots(figsize=(12, 9), dpi=100)
        self.lienzo = FigureCanvasTkAgg(self.figura, master=marco_plot)
        self.lienzo.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        barra_herramientas = NavigationToolbar2Tk(self.lienzo, marco_plot)
        barra_herramientas.update()

        self.actualizar_grafico()

    def cerrar_aplicacion(self):
        self.raiz.quit()
        self.raiz.destroy()
        sys.exit()

    def crear_input(self, padre, texto, valor):
        marco = ttk.Frame(padre)
        marco.pack(fill=tk.X, pady=2)
        ttk.Label(marco, text=texto).pack(anchor=tk.W)
        ent = ttk.Entry(marco)
        ent.insert(0, str(valor))
        ent.pack(fill=tk.X)
        return ent

    def obtener_params(self):
        try:
            return (float(self.ent_alt.get()), float(self.ent_tmin.get()),
                    float(self.ent_tmax.get()), float(self.ent_wmax.get()))
        except ValueError:
            return None

    def abrir_calculadora(self):
        params = self.obtener_params()
        if not params:
            return
        P_MAR = 101325
        PRESION = P_MAR * np.exp(-0.0001255 * params[0])
        CalculadoraPunto(self.raiz, PRESION, self.dibujar_punto)

    def dibujar_punto(self, t, w):
        self.punto_actual = (t, w)
        self.actualizar_grafico()

    def guardar_grafico(self):
        opcion = self.combo_tamano.get()
        if not opcion in OPCIONES_GUARDADO:
            return
        ancho, alto, dpi_g = OPCIONES_GUARDADO[opcion]
        ruta = filedialog.asksaveasfilename(defaultextension=".png", filetypes=[
                                            ("PNG", "*.png"), ("PDF", "*.pdf")])
        if ruta:
            tam_orig = self.figura.get_size_inches()
            try:
                self.figura.set_size_inches(ancho, alto)
                self.figura.savefig(ruta, dpi=dpi_g, bbox_inches='tight')
            finally:
                self.figura.set_size_inches(tam_orig)
                self.lienzo.draw()

    def actualizar_grafico(self):
        params = self.obtener_params()
        if not params:
            return
        ALTITUD, T_MIN, T_MAX, W_MAX_VAL = params
        P_MAR = 101325
        PRESION = P_MAR * np.exp(-0.0001255 * ALTITUD)

        self.ejes.clear()
        T_VEC = np.linspace(T_MIN - 10, T_MAX + 10, 1000)

        # 1. GEOMETRIA ESCALA ENTALPIA (MARCO)
        W_TOP_KG = W_MAX_VAL / 1000.0
        t_corner = obtener_t_saturacion_en_w(W_TOP_KG, PRESION)
        p_corner = np.array([t_corner, W_MAX_VAL])
        w_inicio = razon_humedad_saturacion(T_MIN, PRESION) * 1000.0
        p_inicio = np.array([T_MIN, w_inicio])

        dy = p_corner[1] - p_inicio[1]
        dx = p_corner[0] - p_inicio[0]
        if dx == 0:
            dx = 0.001
        m_diag = dy / dx
        b_diag = p_inicio[1] - m_diag * p_inicio[0]

        vec_diag = np.array([1, m_diag]) / np.linalg.norm([1, m_diag])
        vec_normal = np.array([-vec_diag[1], vec_diag[0]])

        # Dibujo Marco
        self.ejes.plot([p_inicio[0], p_corner[0]], [
                       p_inicio[1], p_corner[1]], color=COLOR_ENTALPIA, linewidth=2.0, zorder=15)
        if t_corner < T_MAX:
            self.ejes.plot([t_corner, T_MAX], [W_MAX_VAL, W_MAX_VAL],
                           color=COLOR_ENTALPIA, linewidth=2.0, zorder=15)

        self.ejes.text(max(t_corner, T_MIN) + 1, W_MAX_VAL + 5, "Entalpia (kJ/kg)",
                       color=COLOR_ENTALPIA, fontweight='bold', fontsize=11, ha='left')

        # Curva Saturacion
        t_sat_plot = np.linspace(T_MIN, t_corner, 300)
        w_sat_plot = razon_humedad_saturacion(t_sat_plot, PRESION) * 1000.0
        self.ejes.plot(t_sat_plot[t_sat_plot >= T_MIN], w_sat_plot[t_sat_plot >= T_MIN],
                       color=COLOR_SATURACION, linewidth=2.5, zorder=20)

        # 2. LINEAS ENTALPIA (Desde Marco hacia Adentro)
        for h in np.arange(-20, 400, PASO_ENTALPIA):
            # a. Calcular interseccion con el marco (Origen del tick)
            t_hit_top = resolver_t_para_h_y_w(h, W_TOP_KG)
            start_point = None
            tick_type = ""

            if t_hit_top >= t_corner:  # Golpea arriba
                if t_hit_top <= T_MAX + 5:
                    start_point = np.array([t_hit_top, W_MAX_VAL])
                    tick_type = 'top'
            else:  # Golpea diagonal
                # A*T^2 + B*T + C = 0 para interseccion recta diagonal vs isentalpica
                m, b = m_diag/1000.0, b_diag/1000.0
                A = 1.86 * m
                B = 2501.0 * m + 1.86 * b + 1.006
                C = 2501.0 * b - h
                raices = np.roots([A, B, C])
                for r in raices:
                    if np.isreal(r) and T_MIN-20 < r < t_corner+5:
                        t_hit = np.real(r)
                        start_point = np.array([t_hit, m_diag*t_hit + b_diag])
                        tick_type = 'diag'
                        break

            if start_point is None:
                continue

            # b. Calcular destino en el grafico (para trazar la linea hacia adentro)
            # Se extiende hasta T_MAX o W=0
            t_line = np.linspace(start_point[0], T_MAX + 30, 100)
            w_line = (h - 1.006 * t_line) / (2501.0 + 1.86 * t_line) * 1000.0

            # c. Filtrar (solo dibujar dentro del grafico y bajo el marco)
            w_diag_lim = m_diag * t_line + b_diag
            mask = (w_line <= W_MAX_VAL+0.05) & (w_line <= w_diag_lim+0.05) & \
                   (w_line >= 0) & (t_line >= T_MIN) & (t_line <= T_MAX)

            if np.any(mask):
                self.ejes.plot(t_line[mask], w_line[mask], color=COLOR_ENTALPIA,
                               linewidth=0.6, alpha=0.8, linestyle='--')

            # d. Dibujar Tick en el marco
            len_tick = 1.2
            if tick_type == 'top':
                p_end = np.array([start_point[0], start_point[1] + len_tick])
                self.ejes.plot([start_point[0], p_end[0]], [start_point[1], p_end[1]],
                               color=COLOR_ENTALPIA, linewidth=1.0, zorder=21)
                if h % PASO_ENTALPIA == 0:
                    fw = 'bold' if h % (PASO_ENTALPIA*2) == 0 else 'normal'
                    if p_end[0] < T_MAX - 1:
                        self.ejes.text(p_end[0], p_end[1] + 0.5, f"{int(h)}",
                                       color=COLOR_ENTALPIA, ha='center', va='bottom', fontsize=7, fontweight=fw, zorder=22)
            elif tick_type == 'diag':
                p_end = start_point + vec_normal * len_tick
                # Solo dibujar si el inicio esta visible en X
                if start_point[0] >= T_MIN:
                    self.ejes.plot([start_point[0], p_end[0]], [start_point[1], p_end[1]],
                                   color=COLOR_ENTALPIA, linewidth=1.0, zorder=21)
                    if h % PASO_ENTALPIA == 0:
                        fw = 'bold' if h % (PASO_ENTALPIA*2) == 0 else 'normal'
                        pos_txt = p_end + vec_normal * 0.5
                        if pos_txt[0] > T_MIN + 1 and pos_txt[1] < W_MAX_VAL + 10:
                            self.ejes.text(pos_txt[0], pos_txt[1], f"{int(h)}",
                                           color=COLOR_ENTALPIA, ha='center', va='center', fontsize=7, fontweight=fw, zorder=22)

        # 3. OTRAS CURVAS
        t_lbl_hr = T_MIN + (T_MAX - T_MIN)*0.5
        w_sat_full = razon_humedad_saturacion(T_VEC, PRESION) * 1000
        for rh in [10, 20, 30, 40, 50, 60, 70, 80, 90]:
            w_rh = razon_humedad_desde_hr(T_VEC, rh/100.0, PRESION) * 1000.0
            w_rh = np.minimum(w_rh, w_sat_full)
            mask = (w_rh <= W_MAX_VAL) & (T_VEC <= T_MAX) & (
                w_rh >= 0) & (T_VEC >= T_MIN)
            self.ejes.plot(T_VEC[mask], w_rh[mask],
                           color=COLOR_HR, linewidth=0.8, alpha=0.6)
            if np.any(mask):
                idx = np.argmin(np.abs(T_VEC[mask] - t_lbl_hr))
                self.ejes.text(T_VEC[mask][idx], w_rh[mask][idx], f"{rh}%",
                               color=COLOR_HR, fontsize=7, va='bottom', ha='center',
                               bbox=dict(facecolor='white', edgecolor='none', pad=0.5, alpha=0.7))

        for twb in np.arange(-10, 60, PASO_TBH):
            if twb < T_MIN or twb > T_MAX:
                continue
            w_s = razon_humedad_saturacion(twb, PRESION)
            h_s = entalpia_aire(twb, w_s)
            w_twb = (h_s - 1.006 * T_VEC) / (2501.0 + 1.86 * T_VEC) * 1000.0
            mask = (w_twb >= 0) & (w_twb <= w_sat_full +
                                   0.1) & (w_twb <= W_MAX_VAL) & (T_VEC >= twb)
            self.ejes.plot(T_VEC[mask], w_twb[mask],
                           color=COLOR_TBH, linestyle=':', linewidth=0.7)
            if len(T_VEC[mask]) > 0:
                self.ejes.text(T_VEC[mask][0], w_twb[mask][0], f"{twb}°C",
                               color=COLOR_TBH, fontsize=6, rotation=-45, ha='right', va='bottom',
                               bbox=dict(facecolor='white', edgecolor='none', pad=0.5, alpha=0.7))

        for v in np.arange(0.70, 1.05, PASO_VOLUMEN):
            w_v = (((v * PRESION)/(287.055 * (T_VEC + 273.15))) - 1) / \
                1.6078 * 1000.0
            mask = (w_v >= 0) & (w_v <= w_sat_full) & (
                w_v <= W_MAX_VAL) & (T_VEC >= T_MIN)
            indices = np.where(mask)[0]
            if len(indices) > 15:
                self.ejes.plot(T_VEC[mask], w_v[mask],
                               color=COLOR_VOL, linewidth=0.5)
                mid = indices[len(indices)//2]
                t_lbl = T_VEC[mid]
                w_lbl = w_v[mid]
                if (t_lbl < T_MAX - 2) and (t_lbl > T_MIN + 2) and (w_lbl < W_MAX_VAL - 1):
                    self.ejes.text(t_lbl, w_lbl, f"{v:.2f} m³/kg", color=COLOR_VOL, fontsize=7, rotation=-65, ha='center',
                                   bbox=dict(facecolor='white', edgecolor='none', pad=0.5, alpha=0.7))

        # 4. PUNTO
        if self.punto_actual:
            pt_t, pt_w = self.punto_actual
            if T_MIN <= pt_t <= T_MAX and 0 <= pt_w <= W_MAX_VAL:
                self.ejes.plot(pt_t, pt_w, 'o', color=COLOR_PUNTO,
                               markersize=6, zorder=30)
                self.ejes.plot([pt_t, pt_t], [0, pt_w],
                               color=COLOR_PUNTO, linestyle=':', linewidth=1)
                self.ejes.plot([pt_t, T_MAX], [pt_w, pt_w],
                               color=COLOR_PUNTO, linestyle=':', linewidth=1)
                self.ejes.text(pt_t + 0.5, pt_w + 0.5, f"({pt_t:.1f}, {pt_w:.1f})",
                               color=COLOR_PUNTO, fontweight='bold', fontsize=8, zorder=31)

        # 5. EJES
        self.ejes.grid(False)
        for tick_w in np.arange(0, W_MAX_VAL + 0.1, 5):
            if tick_w <= 0 or tick_w >= W_MAX_VAL:
                continue
            t_start = max(T_MIN, obtener_t_saturacion_en_w(
                tick_w/1000.0, PRESION))
            self.ejes.plot([t_start, T_MAX], [
                           tick_w, tick_w], color=COLOR_GRILLA, linestyle=':', linewidth=0.5, alpha=0.3)

        x_ticks = np.arange(int(T_MIN/5)*5, T_MAX+1, 5)
        self.ejes.set_xticks(x_ticks)
        for tick_t in x_ticks:
            if tick_t <= T_MIN or tick_t >= T_MAX:
                continue
            w_top = W_MAX_VAL
            if tick_t < t_corner:
                w_top = razon_humedad_saturacion(tick_t, PRESION) * 1000.0
            self.ejes.plot([tick_t, tick_t], [
                           0, w_top], color=COLOR_GRILLA, linestyle=':', linewidth=0.5, alpha=0.3)

        if TEXTO_PIE_PAGINA:
            self.ejes.text(0.99, -0.12, TEXTO_PIE_PAGINA, transform=self.ejes.transAxes,
                           ha='right', va='top', fontsize=8, color='gray')

        self.ejes.set_xlim(T_MIN, T_MAX)
        self.ejes.set_ylim(0, W_MAX_VAL + 10)
        self.ejes.yaxis.tick_right()
        self.ejes.yaxis.set_label_position("right")
        self.ejes.spines['top'].set_visible(False)
        self.ejes.spines['left'].set_visible(False)
        self.ejes.spines['right'].set_bounds(0, W_MAX_VAL)
        self.ejes.set_yticks(np.arange(0, W_MAX_VAL + 1, 5))

        self.ejes.set_xlabel("Temperatura Bulbo Seco (°C)", fontweight='bold')
        self.ejes.set_ylabel("Relacion de Humedad (g/kg)", fontweight='bold')
        self.ejes.set_title(
            f"Carta Psicrometrica ({int(ALTITUD)} msnm)", fontsize=11)
        self.figura.subplots_adjust(bottom=0.15)

        self.lienzo.draw()


raiz = tk.Tk()
raiz.geometry("1200x850")
app = AplicacionPsicrometrica(raiz)
raiz.mainloop()
