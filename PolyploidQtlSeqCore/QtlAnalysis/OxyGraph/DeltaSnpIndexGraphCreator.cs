using OxyPlot;

namespace PolyploidQtlSeqCore.QtlAnalysis.OxyGraph
{
    /// <summary>
    /// ΔSNP-indexグラフクリエーターを作成する。
    /// </summary>
    internal class DeltaSnpIndexGraphCreator : GraphCreator
    {
        private const string TITLE = "ΔSNP-index";
        private const string Y_AXIS_TITLE = "ΔSNP-index";
        private const string FILE_BASENAME = "DeltaSNP-index";

        /// <summary>
        /// ΔSNP-indexグラフクリエーターを作成する。
        /// </summary>
        /// <param name="config">設定</param>
        public DeltaSnpIndexGraphCreator(GraphConfig config)
            : base(config, TITLE, Y_AXIS_TITLE, FILE_BASENAME)
        {
        }

        protected override void AddAnnotations(PlotModel plotModel)
        {
            // アノテーションなし
        }

        protected override void AddSeries(PlotModel plotModel, string chrName, GraphData data)
        {
            var deltaSnpIndexScatterSeries = data.CreateDeltaSnpIndexScatterSeries(chrName);
            plotModel.Series.Add(deltaSnpIndexScatterSeries);

            var p95ThresholdLineSeries = data.CreateP95ThresholdLineSeries(chrName);
            foreach (var series in p95ThresholdLineSeries)
            {
                plotModel.Series.Add(series);
            }

            var p99ThresholdLineSeries = data.CreateP99ThresholdLineSeries(chrName);
            foreach (var series in p99ThresholdLineSeries)
            {
                plotModel.Series.Add(series);
            }

            var averageDeltaSnpIndexLineSeries = data.CreateAverageDeltaSnpIndexLineSeries(chrName);
            foreach (var series in averageDeltaSnpIndexLineSeries)
            {
                plotModel.Series.Add(series);
            }
        }
    }
}
