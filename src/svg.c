/*
    Copyright (C) 2015 Tomas Flouri

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as
    published by the Free Software Foundation, either version 3 of the
    License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.

    You should have received a copy of the GNU Affero General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

    Contact: Tomas Flouri <Tomas.Flouri@h-its.org>,
    Heidelberg Institute for Theoretical Studies,
    Schloss-Wolfsbrunnenweg 35, D-69118 Heidelberg, Germany
*/

#include "delimit.h"

static double scaler = 0;

static long legend_spacing = 10;

static FILE * svg_fp;

static double max_font_len = 0;
static double max_tree_len = 0;

static double canvas_width;

static const char * speciation_color = "#31a354";
static const char * coalesence_color = "#ff0000";

static const char * current_color;

typedef struct coord_s 
{
  double x;
  double y;
} coord_t; 

static coord_t * create_coord(double x, double y)
{
  coord_t * coord = (coord_t *)xmalloc(sizeof(coord_t));
  coord->x = x;
  coord->y = y;
  return coord;
}

static void svg_line(double x1,
                     double y1,
                     double x2,
                     double y2,
                     const char * color,
                     double stroke_width)
{
  fprintf(svg_fp,
          "<line x1=\"%f\" y1=\"%f\" x2=\"%f\" y2=\"%f\" "
          "stroke=\"%s\" stroke-width=\"%f\" />\n",
          x1, y1, x2, y2, color, stroke_width);
}

static void svg_circle(double cx, double cy, double r, const char * color)
{
  fprintf(svg_fp,
          "<circle cx=\"%f\" cy=\"%f\" r=\"%f\" fill=\"%s\" "
          "stroke=\"%s\" />\n",
          cx, cy, r, color, color);
}

static void rtree_set_xcoord(rtree_t * node)
{
  /* create the coordinate info of the node's scaled branch length (edge
     towards root) */
  coord_t * coord = create_coord(node->length * scaler, 0);
  node->data = (void *)coord;

  /* if the node has a parent then add the x coord of the parent such that
     the branch is shifted towards right, otherwise, if the node is the root,
     align it with the left margin */
  if (node->parent)
    coord->x += ((coord_t *)(node->parent->data))->x;
  else
  {
    coord->x = opt_svg_marginleft;
  }

  if (!node->left) 
    return;
  
  /* recursively set coordinates of the other nodes in a pre-order fashion */
  rtree_set_xcoord(node->left);
  rtree_set_xcoord(node->right);
}

static void svg_rtree_plot(rtree_t * node)
{
  double y;
  static int tip_occ = 0;
  double stroke_width = 3;

  /* traverse tree in post-order */
  if (node->left)
  {
    svg_rtree_plot(node->left);
    svg_rtree_plot(node->right);
  }

  if (node->parent)
  {
    double x,px;

    x = ((coord_t *)(node->data))->x;
    px = ((coord_t *)(node->parent->data))->x;

    current_color = (node->parent->event == EVENT_COALESCENT) ? 
                    coalesence_color : speciation_color;

    if (!node->left)
    {
      y = tip_occ * opt_svg_tipspace + opt_svg_margintop + legend_spacing;
      tip_occ++;
    }
    else
    {
      double ly,ry;
      ly = ((coord_t *)(node->left->data))->y;
      ry = ((coord_t *)(node->right->data))->y;
      y = (ly + ry) / 2.0;


      svg_line(x,
               ly,
               x,
               ry,
               (node->event == EVENT_COALESCENT) ?
                 coalesence_color : speciation_color,
               stroke_width);
      
      svg_circle(x,
                 y,
                 opt_svg_inner_radius,
                 (node->event == EVENT_COALESCENT) ?
                   coalesence_color : speciation_color);
    }
    /* horizontal line */
    svg_line(px,
             y,
             x,
             y,
             current_color,
             stroke_width);
    ((coord_t *)(node->data))->y = y;

    if (!node->left)
    {
      fprintf(svg_fp, "<text x=\"%f\" y=\"%f\" "
                      "font-size=\"%ld\" font-family=\"Arial;\">%s</text>\n",
              x+5,
              y+opt_svg_fontsize/3.0,
              opt_svg_fontsize,
              node->label);
    }
    else
      fprintf(svg_fp, "\n");
  }
  else
  {
    double ly,ry,x;
    //    lx = ((coord_t *)(node->left->data))->x;
    ly = ((coord_t *)(node->left->data))->y;
    //    rx = ((coord_t *)(node->right->data))->x;
    ry = ((coord_t *)(node->right->data))->y;
    y = (ly + ry) / 2.0;
    x = opt_svg_marginleft;

    current_color = node->event ? 
                    coalesence_color : speciation_color;

    svg_line(x,
             ly,
             x,
             ry,
             current_color,
             stroke_width);
    svg_circle(x,y,opt_svg_inner_radius, current_color);
  }
}

static void rtree_scaler_init(rtree_t * root)
{
  double len = 0;
  double label_len;
  int i;

  rtree_t ** node_list = (rtree_t **)malloc((size_t)(2 * root->leaves - 1) *
                                             sizeof(rtree_t *));

  rtree_query_tipnodes(root, node_list);

  /* find longest path to root */

  for (i = 0; i < root->leaves; ++i)
  {
    rtree_t * node = node_list[i];

    len = 0;
    while(node)
    {
      len += node->length;
      node = node->parent;
    }
    /* subtract root length */
    len -= root->length;

    if (len > max_tree_len) 
      max_tree_len = len;

    label_len = (opt_svg_fontsize / 1.5) * 
                (node_list[i]->label ? strlen(node_list[i]->label) : 0);

    len = (canvas_width - label_len) / len;
    if (i == 0)
    {
      scaler = len;
      max_font_len = label_len;
    }
    else
      if (len < scaler)
      {
        scaler = len;
        max_font_len = label_len;
      }
  }
  free(node_list);
}

static void svg_rtree_init(rtree_t * root)
{
  long svg_height;

  canvas_width = opt_svg_width - opt_svg_marginleft - opt_svg_marginright;

  /* initialize pixel scaler (scaler) and compute max tree 
     length (max_tree_len) */
  rtree_scaler_init(root);

  svg_height = opt_svg_margintop + legend_spacing + opt_svg_marginbottom + 
               opt_svg_tipspace * root->leaves;


  /* print svg header tag with dimensions and grey border */
  fprintf(svg_fp, "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"%ld\" "
          "height=\"%ld\" style=\"border: 1px solid #cccccc;\">\n",
          opt_svg_width,
          svg_height);

  /* draw legend */
  if (opt_svg_showlegend)
  {
    svg_line(opt_svg_marginleft,
             10,
             (canvas_width - max_font_len)*opt_svg_legend_ratio + 
                                           opt_svg_marginleft,
             10,
             speciation_color,
             3);

    fprintf(svg_fp, "<text x=\"%f\" y=\"%f\" font-size=\"%ld\" "
            "font-family=\"Arial;\">%.*f</text>\n",
            (canvas_width - max_font_len)*opt_svg_legend_ratio + opt_svg_marginleft + 5,
            20-opt_svg_fontsize/3.0,
            (long)opt_svg_fontsize, opt_precision, max_tree_len * opt_svg_legend_ratio);
  }

  /* uncomment to print a dashed border to indicate margins */
  
  /*
  fprintf(svg_fp, "<rect x=\"%ld\" y=\"%ld\" width=\"%ld\" fill=\"none\" "
          "height=\"%ld\" stroke=\"#999999\" stroke-dasharray=\"5,5\" "
          "stroke-width=\"1\" />\n",
          opt_svg_marginleft, 
          opt_svg_margintop + legend_spacing, 
          svg_width - opt_svg_marginleft - opt_svg_marginright,
          svg_height - opt_svg_margintop - legend_spacing - opt_svg_marginbottom);
  */
  
  rtree_set_xcoord(root);

  svg_rtree_plot(root);

  fprintf(svg_fp, "</svg>\n");
}


void cmd_svg(rtree_t * root)
{
  if (!opt_quiet)
    fprintf(stdout,"Creating SVG file %s.txt ...\n", opt_outfile);

  svg_fp = open_file_ext(".svg");

  svg_rtree_init(root);

  fclose(svg_fp);
}
