using OxyPlot;
using OxyPlot.Annotations;

namespace PolyploidQtlSeqCore.QtlAnalysis.OxyGraph
{
    /// <summary>
    /// Windowスコアグラフクリエーター
    /// </summary>
    internal class WindowScoreGraphCreator : GraphCreator
    {
        private const string TITLE = "Window -log10(P)";
        private const string Y_AXIS_TITLE = "-log10(P)";
        private const string FILE_BASENAME = "Window -log10(P)";

        private const int LINE_WIDTH = 2;

        /// <summary>
        /// Windowスコアグラフクリエーターを作成する。
        /// </summary>
        /// <param name="config">設定</param>
        public WindowScoreGraphCreator(GraphConfig config)
            : base(config, TITLE, Y_AXIS_TITLE, FILE_BASENAME)
        {

        }

        protected override void AddAnnotations(PlotModel plotModel)
        {
            var p95LineAnnotation = CreateLineAnnotation(0.05, ColorPalette.P95Color);
            plotModel.Annotations.Add(p95LineAnnotation);

            var p99LineAnnotation = CreateLineAnnotation(0.01, ColorPalette.P99Color);
            plotModel.Annotations.Add(p99LineAnnotation);
        }

        protected override void AddSeries(PlotModel plotModel, string chrName, GraphData data)
        {
            var scoreScatterSeries = data.CreateScoreScatterSeries(chrName);
            plotModel.Series.Add(scoreScatterSeries);
        }

        private static LineAnnotation CreateLineAnnotation(double pValue, OxyColor color)
        {
            var thresholdPValue = new PValue(pValue);
            return new LineAnnotation()
            {
                Type = LineAnnotationType.Horizontal,
                Y = thresholdPValue.ToScore().Value,
                Color = color,
                StrokeThickness = LINE_WIDTH,
                LineStyle = LineStyle.Solid,
            };
        }
    }
}
