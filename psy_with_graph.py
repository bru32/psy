# Simplified chart = PsychrometricChart
# with 
#   graph = NodeGraphEditor
#   and controller = PsychrometricController
#
# Bruce Wernick
# 23 November 2026

import math
import os
import tkinter as tk

if os.name == "nt":
  import ctypes

  try:
    ctypes.windll.shcore.SetProcessDpiAwareness(2)
  except Exception:
    try:
      ctypes.windll.user32.SetProcessDPIAware()
    except Exception:
      pass


# =====================================================================
# 1. PURE PSYCHROMETRIC CHART & EQUATIONS (No Graph Logic)
# =====================================================================
class PsychrometricChart:
  """Draws Tdb vs W axes + saturation curve and provides coordinate mapping."""

  def __init__(self, canvas, width=900, height=650, patm=101.325):
    self.canvas = canvas
    self.width = width
    self.height = height
    self.patm = patm  # Standard atmospheric pressure (kPa)

    # Chart bounds (SI units: °C and kg_w / kg_da)
    self.t_min, self.t_max = 0.0, 50.0
    self.w_min, self.w_max = 0.0, 0.030

    # Canvas plot margins
    self.left = 75
    self.right = width - 40
    self.top = 40
    self.bottom = height - 60

  # --- Moist Air Equations ---
  def saturation_pressure(self, t_db):
    """Magnus-Tetens approximation for saturation vapor pressure (kPa)."""
    return 0.61078 * math.exp((17.27 * t_db) / (t_db + 237.3))

  def saturation_humidity_ratio(self, t_db):
    """Saturation humidity ratio W_sat (kg/kg) at dry bulb temp Tdb (°C)."""
    p_ws = self.saturation_pressure(t_db)
    p_ws = min(p_ws, self.patm * 0.99)
    return 0.621945 * p_ws / (self.patm - p_ws)

  def enthalpy(self, t_db, w):
    """Moist air specific enthalpy h (kJ/kg_da)."""
    return 1.006 * t_db + w * (2501.0 + 1.86 * t_db)

  # --- Coordinate Conversions ---
  def state_to_xy(self, t_db, w):
    """Converts (Tdb, W) to canvas pixel coordinates (x, y)."""
    x_ratio = (t_db - self.t_min) / (self.t_max - self.t_min)
    y_ratio = (w - self.w_min) / (self.w_max - self.w_min)
    x = self.left + x_ratio * (self.right - self.left)
    y = self.bottom - y_ratio * (self.bottom - self.top)
    return x, y

  def xy_to_state(self, x, y):
    """Converts canvas pixel coordinates (x, y) to (Tdb, W)."""
    x_ratio = (x - self.left) / (self.right - self.left)
    y_ratio = (self.bottom - y) / (self.bottom - self.top)
    t_db = self.t_min + x_ratio * (self.t_max - self.t_min)
    w = self.w_min + y_ratio * (self.w_max - self.w_min)
    return t_db, w

  def snap_to_saturation(self, x, y):
    """Constraint hook: snaps (x, y) onto the 100% RH saturation curve."""
    t_db, _ = self.xy_to_state(x, y)
    # Clamp dry bulb so ADP stays within visible chart saturation bounds
    t_db = max(self.t_min, min(31.8, t_db))
    w_sat = self.saturation_humidity_ratio(t_db)
    return self.state_to_xy(t_db, w_sat)

  # --- Background Rendering ---
  def draw(self):
    """Renders axes, gridlines, and the 100% RH saturation curve."""
    # Plot border
    self.canvas.create_rectangle(
        self.left, self.top, self.right, self.bottom, outline="#adb5bd", width=1
    )

    # X-Axis: Dry Bulb Temperature (°C)
    for t in range(int(self.t_min), int(self.t_max) + 1, 5):
      x, _ = self.state_to_xy(t, self.w_min)
      self.canvas.create_line(
          x, self.top, x, self.bottom, fill="#e9ecef", dash=(2, 2)
      )
      self.canvas.create_line(x, self.bottom, x, self.bottom + 5, fill="#495057")
      self.canvas.create_text(
          x, self.bottom + 18, text=f"{t}°C", font=("Arial", 8), fill="#495057"
      )

    self.canvas.create_text(
        (self.left + self.right) / 2,
        self.bottom + 40,
        text="Dry Bulb Temperature (°C)",
        font=("Arial", 10, "bold"),
        fill="#212529",
    )

    # Y-Axis: Humidity Ratio (kg/kg)
    steps = 6
    for i in range(steps + 1):
      w = self.w_min + i * (self.w_max - self.w_min) / steps
      _, y = self.state_to_xy(self.t_min, w)
      self.canvas.create_line(
          self.left, y, self.right, y, fill="#e9ecef", dash=(2, 2)
      )
      self.canvas.create_line(self.left - 5, y, self.left, y, fill="#495057")
      self.canvas.create_text(
          self.left - 35,
          y,
          text=f"{w:.3f}",
          font=("Arial", 8),
          fill="#495057",
      )

    self.canvas.create_text(
        20,
        (self.top + self.bottom) / 2,
        text="Humidity Ratio W (kg/kg)",
        angle=90,
        font=("Arial", 10, "bold"),
        fill="#212529",
    )

    # 100% RH Saturation Curve
    curve_pts = []
    for i in range(201):
      t = self.t_min + (i / 200.0) * (self.t_max - self.t_min)
      w_sat = self.saturation_humidity_ratio(t)
      if w_sat <= self.w_max:
        px, py = self.state_to_xy(t, w_sat)
        curve_pts.extend([px, py])

    if len(curve_pts) >= 4:
      self.canvas.create_line(
          curve_pts, fill="#0077b6", width=2.5, smooth=True, tags="chart_bg"
      )


# =====================================================================
# 2. PURE ABSTRACT GRAPH EDITOR (No Psychrometric Knowledge)
# =====================================================================
class NodeGraphEditor:
  """Manages nodes, MIX links, ADP links, and optional position constraints."""

  def __init__(self, canvas, root):
    self.canvas = canvas
    self.root = root

    self.radius = 12
    self.slider_radius = 9

    # item_id -> {
    #   "x": float, "y": float,
    #   "parent_link": link_dict or None,
    #   "role": "STANDARD" | "MIX" | "ADP",
    #   "constraint": callable(x, y) -> (x, y) or None
    # }
    self.nodes = {}
    self.links = []
    self.selected_node = None
    self.drag_data = {"type": None, "item": None, "x": 0, "y": 0}

    # External Hooks for Controller
    self.on_node_created = None  # callback(node_id, role)
    self.on_graph_changed = None  # callback()

    # Bindings
    self.canvas.bind("<ButtonPress-1>", self.on_left_down)
    self.canvas.bind("<B1-Motion>", self.on_left_drag)
    self.canvas.bind("<ButtonRelease-1>", self.on_left_up)
    self.canvas.bind("<ButtonPress-3>", self.on_right_click)
    self.canvas.bind("<Shift-ButtonPress-3>", self.on_shift_right_click)
    self.root.bind("<Delete>", self.on_delete)
    self.root.bind("<BackSpace>", self.on_delete)

  def notify_changed(self):
    if self.on_graph_changed:
      self.on_graph_changed()

  def get_item_at(self, x, y, tag_filter):
    overlapping = self.canvas.find_overlapping(x - 4, y - 4, x + 4, y + 4)
    for item in reversed(overlapping):
      if tag_filter in self.canvas.gettags(item):
        return item
    return None

  def select_node(self, node_id):
    if self.selected_node and self.selected_node in self.nodes:
      self.canvas.itemconfig(self.selected_node, outline="black", width=2)
    self.selected_node = node_id
    if self.selected_node and self.selected_node in self.nodes:
      self.canvas.itemconfig(self.selected_node, outline="#0066cc", width=3)

  def create_node(self, x, y, parent_link=None, role="STANDARD", fill="white"):
    r = self.slider_radius if parent_link else self.radius
    node_id = self.canvas.create_oval(
        x - r,
        y - r,
        x + r,
        y + r,
        fill=fill,
        outline="black",
        width=2,
        tags="node",
    )
    self.nodes[node_id] = {
        "x": x,
        "y": y,
        "parent_link": parent_link,
        "role": role,
        "constraint": None,
    }
    self.select_node(node_id)

    # Notify controller so it can attach constraints (e.g., saturation snap)
    if self.on_node_created:
      self.on_node_created(node_id, role)

    # Apply constraint immediately if one was attached
    self.set_node_position(node_id, x, y)
    self.notify_changed()
    return node_id

  def set_node_position(self, node_id, x, y):
    """Applies any registered position constraint before placing the node."""
    node = self.nodes[node_id]
    if node["constraint"]:
      x, y = node["constraint"](x, y)

    node["x"] = x
    node["y"] = y
    r = self.slider_radius if node["parent_link"] else self.radius
    self.canvas.coords(node_id, x - r, y - r, x + r, y + r)

  def create_link(self, n1, n2):
    if n1 == n2:
      return
    p1 = self.nodes[n1]["parent_link"]
    p2 = self.nodes[n2]["parent_link"]
    if (p1 and n2 in (p1["n1"], p1["n2"])) or (
        p2 and n1 in (p2["n1"], p2["n2"])
    ):
      return
    for link in self.links:
      if {link["n1"], link["n2"]} == {n1, n2}:
        return

    x1, y1 = self.nodes[n1]["x"], self.nodes[n1]["y"]
    x2, y2 = self.nodes[n2]["x"], self.nodes[n2]["y"]
    line_id = self.canvas.create_line(
        x1, y1, x2, y2, fill="#333333", width=2, tags="link"
    )
    # Keep links above chart background but below nodes
    self.canvas.tag_raise(line_id, "chart_bg")

    self.links.append({
        "id": line_id,
        "n1": n1,
        "n2": n2,
        "type": "STANDARD",
        "slider_node": None,
        "t": 0.5,
        "bypass_factor": 0.0,
        "label_id": None,
    })
    self.notify_changed()

  def project_raw_t(self, px, py, link):
    x1, y1 = self.nodes[link["n1"]]["x"], self.nodes[link["n1"]]["y"]
    x2, y2 = self.nodes[link["n2"]]["x"], self.nodes[link["n2"]]["y"]
    vx, vy = x2 - x1, y2 - y1
    seg_len_sq = vx * vx + vy * vy
    return (
        0.5
        if seg_len_sq == 0
        else ((px - x1) * vx + (py - y1) * vy) / seg_len_sq
    )

  def project_onto_segment(self, px, py, x1, y1, x2, y2):
    vx, vy = x2 - x1, y2 - y1
    seg_len_sq = vx * vx + vy * vy
    if seg_len_sq == 0:
      return 0.5
    u = ((px - x1) * vx + (py - y1) * vy) / seg_len_sq
    return max(0.05, min(0.95, u))

  def get_adp_link_for_offcoil(self, node_id):
    for link in self.links:
      if link["type"] == "ADP" and link["n2"] == node_id:
        return link
    return None

  def remove_link(self, link):
    if link not in self.links:
      return
    self.links.remove(link)
    slider_id = link["slider_node"]
    if slider_id and slider_id in self.nodes:
      self.delete_node_cascade(slider_id)
    if link["label_id"]:
      self.canvas.delete(link["label_id"])
    self.canvas.delete(link["id"])

  def delete_node_cascade(self, node_id):
    if node_id not in self.nodes:
      return
    if self.selected_node == node_id:
      self.selected_node = None
    for link in [
        l for l in self.links if l["n1"] == node_id or l["n2"] == node_id
    ]:
      self.remove_link(link)
    self.canvas.delete(node_id)
    del self.nodes[node_id]

  def clear_link_slider(self, link):
    slider_id = link["slider_node"]
    link["slider_node"] = None
    if slider_id and slider_id in self.nodes:
      self.delete_node_cascade(slider_id)
    if link["label_id"]:
      self.canvas.delete(link["label_id"])
    link["label_id"] = None
    link["type"] = "STANDARD"
    link["bypass_factor"] = 0.0
    self.canvas.itemconfig(link["id"], fill="#333333", width=2)
    self.update_links()

  def set_link_mode(self, link, mode, click_x=None, click_y=None):
    self.clear_link_slider(link)
    link["type"] = mode

    if mode == "MIX":
      t = (
          self.project_raw_t(click_x, click_y, link)
          if click_x is not None
          else 0.5
      )
      link["t"] = max(0.0, min(1.0, t))
      link["bypass_factor"] = 0.0
      fill_color = "#ffeaa7"
      self.canvas.itemconfig(link["id"], fill="#d35400", width=2.5)

    elif mode == "ADP":
      if click_x is not None and click_y is not None:
        if self.project_raw_t(click_x, click_y, link) < 0.5:
          link["n1"], link["n2"] = link["n2"], link["n1"]
      link["t"] = 1.35
      link["bypass_factor"] = 0.35
      fill_color = "#81ecec"
      self.canvas.itemconfig(link["id"], fill="#008080", width=2.5)

    x1, y1 = self.nodes[link["n1"]]["x"], self.nodes[link["n1"]]["y"]
    x2, y2 = self.nodes[link["n2"]]["x"], self.nodes[link["n2"]]["y"]
    sx = x1 + link["t"] * (x2 - x1)
    sy = y1 + link["t"] * (y2 - y1)

    slider_id = self.create_node(
        sx, sy, parent_link=link, role=mode, fill=fill_color
    )
    link["slider_node"] = slider_id
    link["label_id"] = self.canvas.create_text(
        sx, sy - 18, text="", font=("Arial", 8, "bold"), fill="#2d3436"
    )

    # If ADP was snapped onto a constraint (like the saturation curve), align Off-Coil n2
    if mode == "ADP":
      snapped_x = self.nodes[slider_id]["x"]
      snapped_y = self.nodes[slider_id]["y"]
      new_x2 = x1 + (snapped_x - x1) / link["t"]
      new_y2 = y1 + (snapped_y - y1) / link["t"]
      self.set_node_position(link["n2"], new_x2, new_y2)

    self.update_links()

  def move_node(self, node_id, dx, dy):
    if self.nodes[node_id]["parent_link"] is not None:
      return
    nx = self.nodes[node_id]["x"] + dx
    ny = self.nodes[node_id]["y"] + dy
    self.set_node_position(node_id, nx, ny)

  def update_links(self):
    for _ in range(len(self.links) + 1):
      moved = False
      for link in self.links:
        n1, n2 = link["n1"], link["n2"]
        x1, y1 = self.nodes[n1]["x"], self.nodes[n1]["y"]
        x2, y2 = self.nodes[n2]["x"], self.nodes[n2]["y"]

        slider_id = link["slider_node"]
        if slider_id and slider_id in self.nodes:
          if link["type"] == "MIX":
            t = link["t"]
            sx = x1 + t * (x2 - x1)
            sy = y1 + t * (y2 - y1)
            if (
                abs(sx - self.nodes[slider_id]["x"]) > 1e-4
                or abs(sy - self.nodes[slider_id]["y"]) > 1e-4
            ):
              moved = True
            self.set_node_position(slider_id, sx, sy)
            self.canvas.coords(link["id"], x1, y1, x2, y2)
            label_text = f"mix: {t:.2f}"

          elif link["type"] == "ADP":
            # Keep ADP slider on its constraint and align Off-Coil n2 along n1 -> ADP
            sx = self.nodes[slider_id]["x"]
            sy = self.nodes[slider_id]["y"]
            t = link["t"]
            new_x2 = x1 + (sx - x1) / t
            new_y2 = y1 + (sy - y1) / t
            if (
                abs(new_x2 - x2) > 1e-4
                or abs(new_y2 - y2) > 1e-4
            ):
              moved = True
            self.set_node_position(n2, new_x2, new_y2)
            self.canvas.coords(link["id"], x1, y1, sx, sy)
            label_text = f"ADP (bf: {link['bypass_factor']:.2f})"

          if link["label_id"]:
            self.canvas.coords(link["label_id"], sx, sy - 18)
            self.canvas.itemconfig(link["label_id"], text=label_text)
        else:
          self.canvas.coords(link["id"], x1, y1, x2, y2)

      if not moved:
        break

    self.notify_changed()

  # --- Mouse & Keyboard Handlers ---
  def on_left_down(self, event):
    node = self.get_item_at(event.x, event.y, "node")
    if node:
      self.select_node(node)
      self.drag_data = {
          "type": "node",
          "item": node,
          "x": event.x,
          "y": event.y,
      }
      return

    link_id = self.get_item_at(event.x, event.y, "link")
    if link_id:
      link_obj = next((l for l in self.links if l["id"] == link_id), None)
      if link_obj:
        self.drag_data = {
            "type": "link",
            "item": link_obj,
            "x": event.x,
            "y": event.y,
        }
        return

    self.create_node(event.x, event.y)

  def on_left_drag(self, event):
    if not self.drag_data["type"]:
      return
    dx = event.x - self.drag_data["x"]
    dy = event.y - self.drag_data["y"]
    self.drag_data["x"] = event.x
    self.drag_data["y"] = event.y

    if self.drag_data["type"] == "node":
      node_id = self.drag_data["item"]
      parent_link = self.nodes[node_id]["parent_link"]

      if parent_link:
        if parent_link["type"] == "MIX":
          raw_t = self.project_raw_t(event.x, event.y, parent_link)
          parent_link["t"] = max(0.0, min(1.0, raw_t))
          self.update_links()

        elif parent_link["type"] == "ADP":
          # set_node_position automatically calls the saturation constraint hook!
          self.set_node_position(node_id, float(event.x), float(event.y))
          self.update_links()
      else:
        adp_link = self.get_adp_link_for_offcoil(node_id)
        if adp_link:
          n1 = adp_link["n1"]
          adp_id = adp_link["slider_node"]
          x1, y1 = self.nodes[n1]["x"], self.nodes[n1]["y"]
          sx, sy = self.nodes[adp_id]["x"], self.nodes[adp_id]["y"]

          u = self.project_onto_segment(event.x, event.y, x1, y1, sx, sy)
          t = 1.0 / u
          adp_link["t"] = t
          adp_link["bypass_factor"] = t - 1.0
          self.set_node_position(
              node_id, x1 + u * (sx - x1), y1 + u * (sy - y1)
          )
          self.update_links()
        else:
          self.move_node(node_id, dx, dy)
          self.update_links()

    elif self.drag_data["type"] == "link":
      link_obj = self.drag_data["item"]
      self.move_node(link_obj["n1"], dx, dy)
      self.move_node(link_obj["n2"], dx, dy)
      self.update_links()

    self.canvas.update_idletasks()

  def on_left_up(self, event):
    self.drag_data = {"type": None, "item": None, "x": 0, "y": 0}

  def on_right_click(self, event):
    target_node = self.get_item_at(event.x, event.y, "node")
    if target_node:
      if self.selected_node and self.selected_node in self.nodes:
        self.create_link(self.selected_node, target_node)
      self.select_node(target_node)
      return

    target_link_id = self.get_item_at(event.x, event.y, "link")
    if target_link_id:
      link_obj = next((l for l in self.links if l["id"] == target_link_id), None)
      if link_obj:
        if link_obj["type"] == "STANDARD":
          self.set_link_mode(link_obj, "MIX", event.x, event.y)
        elif link_obj["type"] == "MIX":
          self.set_link_mode(link_obj, "ADP", event.x, event.y)
        else:
          self.clear_link_slider(link_obj)

  def on_shift_right_click(self, event):
    target_link_id = self.get_item_at(event.x, event.y, "link")
    if target_link_id:
      link_obj = next((l for l in self.links if l["id"] == target_link_id), None)
      if link_obj:
        self.set_link_mode(link_obj, "ADP", event.x, event.y)

  def on_delete(self, event):
    if not self.selected_node or self.selected_node not in self.nodes:
      return
    node_to_delete = self.selected_node
    parent_link = self.nodes[node_to_delete]["parent_link"]
    if parent_link:
      self.clear_link_slider(parent_link)
    else:
      self.delete_node_cascade(node_to_delete)
    self.notify_changed()


# =====================================================================
# 3. CONTROLLER (Bridges Graph Events & Psychrometric Math)
# =====================================================================
class PsychrometricController:
  """Connects the abstract graph to the psychrometric chart without coupling them."""

  def __init__(self, chart, graph, status_var):
    self.chart = chart
    self.graph = graph
    self.status_var = status_var

    # Subscribe to graph hooks
    self.graph.on_node_created = self.handle_node_created
    self.graph.on_graph_changed = self.handle_graph_changed

  def handle_node_created(self, node_id, role):
    """Attaches the saturation curve snap constraint to any ADP node."""
    if role == "ADP":
      self.graph.nodes[node_id]["constraint"] = self.chart.snap_to_saturation

  def handle_graph_changed(self):
    """Translates the currently selected node's pixel coordinates into moist air properties."""
    sel = self.graph.selected_node
    if sel and sel in self.graph.nodes:
      x = self.graph.nodes[sel]["x"]
      y = self.graph.nodes[sel]["y"]
      role = self.graph.nodes[sel]["role"]

      t_db, w = self.chart.xy_to_state(x, y)
      h = self.chart.enthalpy(t_db, w)
      self.status_var.set(
          f"Selected [{role}]  |  Tdb: {t_db:.2f} °C  |  "
          f"W: {w:.5f} kg/kg  |  Enthalpy h: {h:.2f} kJ/kg"
      )
    else:
      self.status_var.set(
          "Left-click: Add Node | Right-click Node: Link | "
          "Right-click Link: MIX | Shift+Right-click Link: ADP"
      )


  def seed_sample_hvac_cycle(self):
    """Places a pre-built Outside Air -> Mix -> Cooling Coil -> Room Air cycle."""
    # 1. Outside Air (35°C, 0.018 kg/kg) & Room Air (24°C, 0.0095 kg/kg)
    oa_x, oa_y = self.chart.state_to_xy(35.0, 0.018)
    ra_x, ra_y = self.chart.state_to_xy(24.0, 0.0095)
    oa_id = self.graph.create_node(oa_x, oa_y)
    ra_id = self.graph.create_node(ra_x, ra_y)

    # 2. Create MIX link between OA and RA (30% Outside Air)
    self.graph.create_link(oa_id, ra_id)
    mix_link = self.graph.links[-1]
    mix_x, mix_y = self.chart.state_to_xy(27.3, 0.01205)
    self.graph.set_link_mode(mix_link, "MIX", mix_x, mix_y)
    mix_node_id = mix_link["slider_node"]

    # 3. Create Off-Coil node & ADP link from MIX node (ADP snaps to saturation curve!)
    oc_x, oc_y = self.chart.state_to_xy(14.0, 0.0090)
    oc_id = self.graph.create_node(oc_x, oc_y)
    self.graph.create_link(mix_node_id, oc_id)
    adp_link = self.graph.links[-1]
    self.graph.set_link_mode(adp_link, "ADP", oc_x, oc_y)

    # 4. Complete the cycle: Off-Coil -> Room Air (Room Heat Load line)
    self.graph.create_link(oc_id, ra_id)
    self.graph.select_node(adp_link["slider_node"])
    self.handle_graph_changed()


if __name__ == "__main__":

  root = tk.Tk()
  root.title("Psychrometric Cycle Designer (Clean MVC Separation)")

  status_var = tk.StringVar()
  status_bar = tk.Label(root,
      textvariable=status_var,
      bd=1, relief=tk.SUNKEN,
      anchor=tk.W, font=("Arial", 10),
      bg="#e9ecef", padx=10, pady=6)
  status_bar.pack(side=tk.BOTTOM, fill=tk.X)

  canvas = tk.Canvas(root, width=900, height=650, bg="white")
  canvas.pack(fill=tk.BOTH, expand=True)

  # Instantiate the 3 decoupled layers
  chart = PsychrometricChart(canvas, width=900, height=650)
  chart.draw()

  graph = NodeGraphEditor(canvas, root)
  controller = PsychrometricController(chart, graph, status_var)

  # Load a sample cycle so you can immediately drag the ADP node along the saturation curve
  controller.seed_sample_hvac_cycle()

  root.mainloop()
