using OxyPlot;
using OxyPlot.Annotations;

namespace PolyploidQtlSeqCore.QtlAnalysis.OxyGraph
{
    /// <summary>
    /// QTL数グラフクリエーター
    /// </summary>
    internal class QtlVariantCountGraphCreator : GraphCreator
    {
        private const string TITLE = "QTL Variant Count";
        private const string Y_AXIS_TITLE = "QTL Count";
        private const string FILE_BASENAME = "QTL Variant Count";

        /// <summary>
        /// QTL数グラフクリエーターを作成する。
        /// </summary>
        /// <param name="config">設定</param>
        public QtlVariantCountGraphCreator(GraphConfig config)
            : base(config, TITLE, Y_AXIS_TITLE, FILE_BASENAME)
        {
        }

        protected override void AddAnnotations(PlotModel plotModel)
        {
            var zeroLineAnnotation = new LineAnnotation()
            {
                Type = LineAnnotationType.Horizontal,
                Y = 0,
                Color = OxyColors.Black,
                StrokeThickness = 2,
                LineStyle = LineStyle.Solid,
            };

            plotModel.Annotations.Add(zeroLineAnnotation);
        }

        protected override void AddSeries(PlotModel plotModel, string chrName, GraphData data)
        {
            var p95QtlCountScatterSeries = data.CreateP95QtlCountScatterSeries(chrName);
            plotModel.Series.Add(p95QtlCountScatterSeries);

            var p99QtlCountScatterSeries = data.CreateP99QtlCountScatterSeries(chrName);
            plotModel.Series.Add(p99QtlCountScatterSeries);
        }
    }
}
