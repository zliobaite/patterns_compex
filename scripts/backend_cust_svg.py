import numpy
import re
import six

# from matplotlib import verbose, __version__, rcParams
from matplotlib import __version__, rcParams

import matplotlib.backend_bases
from matplotlib.backends.backend_mixed import MixedModeRenderer
from matplotlib.backends.backend_svg import *  # RendererPgf, FigureCanvasPgf
from matplotlib.text import Text

import pdb


class RendererSVGCust(RendererSVG):

    # def __init__(self, width, height, svgwriter, basename=None, image_dpi=72, fig=None):
    #     RendererSVG.__init__(self, width, height, svgwriter, basename, image_dpi)
    #     self.fig = fig

    # def set_figure(self, fig):
    #     self.fig
    # def get_figure(self):
    #     return self.fig

    # def open_group(self, s, gid=None):
    #     """
    #     Open a grouping element with label *s*. If *gid* is given, use
    #     *gid* as the id of the group.
    #     """
    #     if gid:
    #         if gid < 0:
    #             self.writer.start('g', attrib={"id": "%s_%d" % (s, -gid), "class": 'lbl_txt'})
    #         else:
    #             self.writer.start('g', id=gid)
    #     else:
    #         self._groupd[s] = self._groupd.get(s, 0) + 1
    #         self.writer.start('g', id="%s_%d" % (s, self._groupd[s]))

    def _write_default_style(self):
        writer = self.writer
        default_style = generate_css({
            'stroke-linejoin': 'round',
            'stroke-linecap': 'butt'})
        writer.start('defs')
        writer.start('style', type='text/css')
        writer.data('*{%s}\n' % default_style)
        writer.data('g.lbl text {%s}\n' % generate_css({'display': 'none', 'font-size': 'small', 'font-style': 'italic', 'font-family': 'Arial, Helvetica, sans-serif'}))
        writer.data('g.lbl:hover text {%s}\n' % generate_css({'display': 'block'}))
        writer.end('style')
        writer.start('filter', {"x": "0", "y": "0", "width": "100%", "height": "100%", "id": "bg-lbl"})
        writer.start('feFlood', {"flood-color": "white"})
        writer.end('feFlood')
        writer.start('feGaussianBlur', {"stdDeviation": "2"})
        writer.end('feGaussianBlur')
        writer.end('filter')
        writer.end('defs')

    def text_lbl(self, writer, lbl, x, y):
        #self.writer.start('span', attrib={"class": 'lbl_span'})
        writer.element('text', lbl,
                       attrib={'x': short_float_fmt(x), 'y': short_float_fmt(y),
                               'filter': "url(#bg-lbl)"})
        writer.element('text', lbl,
                       attrib={'x': short_float_fmt(x), 'y': short_float_fmt(y)})

        # self.writer.end('span')

    # def text_lbl(self, writer, lbl, path, transform):
    #     # trans_and_flip = self._make_flip_transform(transform)
    #     lbl_bb = path.get_extents().get_points()
    #     x = (lbl_bb[0,0]+lbl_bb[1,0])/2
    #     y = (lbl_bb[0,1]+lbl_bb[1,1])/2

    #     print(x, y, lbl)
    #     # writer.element('text', lbl, #
    #     #             attrib={'x': short_float_fmt(x), 'y': short_float_fmt(y)})
    #     textobj = Text(x, y, "XX"+lbl)
    #     textobj.set_figure(self.get_figure())
    #     s = 'text'
    #     self._groupd[s] = self._groupd.get(s, 0) + 1
    #     textobj.set_gid(-self._groupd[s])
    #     textobj.draw(self)
    #     # if rcParams['svg.fonttype'] == 'path':
    #     #     self._draw_text_as_path(gc, x, y, s, prop, angle, ismath, mtext)
    #     # else:
    #     #     self._draw_text_as_text(gc, x, y, s, prop, angle, ismath, mtext)

    # def draw_tex(self, gc, x, y, s, prop, angle, ismath='TeX!', mtext=None):
    #     pdb.set_trace()
    #     self._draw_text_as_path(gc, x, y, s, prop, angle, ismath="TeX")

    def draw_path(self, gc, path, transform, rgbFace=None):
        # pdb.set_trace()
        trans_and_flip = self._make_flip_transform(transform)
        clip = (rgbFace is None and gc.get_hatch_path() is None)
        simplify = path.should_simplify and clip
        path_data = self._convert_path(
            path, trans_and_flip, clip=clip, simplify=simplify,
            sketch=gc.get_sketch_params())

        attrib = {}
        attrib['style'] = self._get_style(gc, rgbFace)

        clipid = self._get_clip(gc)
        if clipid is not None:
            attrib['clip-path'] = 'url(#%s)' % clipid

        lbl = None
        if gc.get_url() is not None:
            if re.match("LBL:", gc.get_url()):
                lbl = gc.get_url()[4:]
                lbl_bb = path.get_extents(trans_and_flip).get_points()
                self.writer.start('g', attrib={'class': 'lbl'})
                self.text_lbl(self.writer, lbl, (lbl_bb[0, 0]+lbl_bb[1, 0])/2, (lbl_bb[0, 1]+lbl_bb[1, 1])/2)

            else:
                self.writer.start('a', {'xlink:href': gc.get_url()})
        self.writer.element('path', d=path_data, attrib=attrib)
        if gc.get_url() is not None:
            if lbl is not None:
                self.writer.end('g')
            else:
                self.writer.end('a')

    def draw_markers(self, gc, marker_path, marker_trans, path, trans, rgbFace=None):
        if not len(path.vertices):
            return

        writer = self.writer
        path_data = self._convert_path(
            marker_path,
            marker_trans + Affine2D().scale(1.0, -1.0),
            simplify=False)
        style = self._get_style_dict(gc, rgbFace)
        dictkey = (path_data, generate_css(style))
        oid = self._markers.get(dictkey)
        style = generate_css({k: v for k, v in six.iteritems(style)
                              if k.startswith('stroke')})

        if oid is None:
            oid = self._make_id('m', dictkey)
            writer.start('defs')
            writer.element('path', id=oid, d=path_data, style=style)
            writer.end('defs')
            self._markers[dictkey] = oid

        attrib = {}
        clipid = self._get_clip(gc)
        url = gc.get_url()
        lbl = None
        if url is not None and re.match("LBL:", url):
            lbl = url[4:]
            attrib['class'] = 'lbl'
        elif clipid is not None:
            attrib['clip-path'] = 'url(#%s)' % clipid

        writer.start('g', attrib=attrib)

        trans_and_flip = self._make_flip_transform(trans)
        if lbl is not None:
            lbl_bb = path.get_extents(trans_and_flip).get_points()
            self.text_lbl(writer, lbl, (lbl_bb[0, 0]+lbl_bb[1, 0])/2, (lbl_bb[0, 1]+lbl_bb[1, 1])/2)

        attrib = {'xlink:href': '#%s' % oid}
        clip = (0, 0, self.width*72, self.height*72)
        for vertices, code in path.iter_segments(
                trans_and_flip, clip=clip, simplify=False):
            if len(vertices):
                x, y = vertices[-2:]
                attrib['x'] = short_float_fmt(x)
                attrib['y'] = short_float_fmt(y)
                attrib['style'] = self._get_style(gc, rgbFace)
                writer.element('use', attrib=attrib)
        writer.end('g')

    def draw_path_collection(self, gc, master_transform, paths, all_transforms,
                             offsets, offsetTrans, facecolors, edgecolors,
                             linewidths, linestyles, antialiaseds, urls,
                             offset_position):
        # Is the optimization worth it? Rough calculation:
        # cost of emitting a path in-line is
        #    (len_path + 5) * uses_per_path
        # cost of definition+use is
        #    (len_path + 3) + 9 * uses_per_path
        # pdb.set_trace()
        len_path = len(paths[0].vertices) if len(paths) > 0 else 0
        uses_per_path = self._iter_collection_uses_per_path(
            paths, all_transforms, offsets, facecolors, edgecolors)
        should_do_optimization = \
            len_path + 9 * uses_per_path + 3 < (len_path + 5) * uses_per_path
        if False:  # not should_do_optimization:
            return RendererBase.draw_path_collection(
                self, gc, master_transform, paths, all_transforms,
                offsets, offsetTrans, facecolors, edgecolors,
                linewidths, linestyles, antialiaseds, urls,
                offset_position)

        writer = self.writer
        path_codes = []
        writer.start('defs')
        for i, (path, transform) in enumerate(self._iter_collection_raw_paths(
                master_transform, paths, all_transforms)):
            transform = Affine2D(transform.get_matrix()).scale(1.0, -1.0)
            d = self._convert_path(path, transform, simplify=False)
            oid = 'C%x_%x_%s' % (self._path_collection_id, i,
                                 self._make_id('', d))
            writer.element('path', id=oid, d=d)
            path_codes.append(oid)
        writer.end('defs')

        for xo, yo, path_id, gc0, rgbFace in self._iter_collection(
                gc, master_transform, all_transforms, path_codes, offsets,
                offsetTrans, facecolors, edgecolors, linewidths, linestyles,
                antialiaseds, urls, offset_position):
            clipid = self._get_clip(gc0)
            url = gc0.get_url()
            lbl = None
            if url is not None:
                if re.match("LBL:", url):
                    lbl = url[4:]
                    attrib = {'class': 'lbl'}
                    # if clipid is not None:
                    #     attrib.update({'clip-path': 'url(#%s)' % clipid})
                    writer.start('g', attrib=attrib)
                    self.text_lbl(writer, lbl, xo, self.height - yo)

                else:
                    writer.start('a', attrib={'xlink:href': url})

            if lbl is None and clipid is not None:
                writer.start('g', attrib={'clip-path': 'url(#%s)' % clipid})
            attrib = {
                'xlink:href': '#%s' % path_id,
                'x': short_float_fmt(xo),
                'y': short_float_fmt(self.height - yo),
                'style': self._get_style(gc0, rgbFace)
            }
            writer.element('use', attrib=attrib)
            if lbl is None and clipid is not None:
                writer.end('g')
            if url is not None:
                if lbl is not None:
                    self.writer.end('g')
                else:
                    writer.end('a')

        self._path_collection_id += 1


class FigureCanvasSVGCust(FigureCanvasSVG):

    def _print_svg(self, filename, svgwriter, **kwargs):
        image_dpi = kwargs.pop("dpi", 72)
        self.figure.set_dpi(72.0)
        width, height = self.figure.get_size_inches()
        w, h = width*72, height*72

        _bbox_inches_restore = kwargs.pop("bbox_inches_restore", None)
        renderer = MixedModeRenderer(
            self.figure,
            width, height, image_dpi, RendererSVGCust(w, h, svgwriter, filename, image_dpi),  # , self.figure),
            bbox_inches_restore=_bbox_inches_restore)

        self.figure.draw(renderer)
        renderer.finalize()


matplotlib.backend_bases.register_backend("svg", FigureCanvasSVGCust, description="custom SVG backend")
