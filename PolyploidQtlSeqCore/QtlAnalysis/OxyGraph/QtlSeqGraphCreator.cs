using PolyploidQtlSeqCore.QtlAnalysis.SlidingWindow;

namespace PolyploidQtlSeqCore.QtlAnalysis.OxyGraph
{
    /// <summary>
    /// QtlSeqグラフクリエーター
    /// </summary>
    internal class QtlSeqGraphCreator
    {
        private readonly GraphSettings _settings;

        /// <summary>
        /// QtlSeqグラフクリエーターを作成する。
        /// </summary>
        /// <param name="setting">グラフ設定</param>
        public QtlSeqGraphCreator(GraphSettings setting)
        {
            _settings = setting;
        }

        /// <summary>
        /// グラフを作成する。
        /// </summary>
        /// <param name="outDir">出力ディレクトリ</param>
        /// <param name="allVariants">全変異データ</param>
        /// <param name="allWindows">全Windowデータ</param>
        public void Create(OutputDirectory outDir, SnpIndexVariantWithSlidingWindowQtl[] allVariants, Window[] allWindows)
        {
            var graphData = new GraphData(allVariants, allWindows, _settings.XAxisMajorStep);
            var graphCreator = GetGrhaphCreators(graphData.GraphAxes);

            foreach (var chrName in graphData.ChrNames)
            {
                var mergeGraphFilePath = outDir.CreateFilePath($"{chrName}.png");

                var files = graphCreator.Select(x => x.Create(outDir, chrName, graphData)).ToArray();
                var graphFiles = new GraphFiles(files);
                graphFiles.VerticalMerge(mergeGraphFilePath, _settings);
                graphFiles.Delete();
            }
        }

        /// <summary>
        /// グラフクリエーター配列を取得する。
        /// </summary>
        /// <param name="graphAxes">グラフ軸</param>
        /// <returns>グラフクリエーター配列</returns>
        private GraphCreator[] GetGrhaphCreators(GraphAxes graphAxes)
        {
            var correctAxies = graphAxes.CorrectYAxis(_settings);

            var snpIndexConfig = correctAxies.CreateSnpIndexGraphConfig(_settings);
            var deltaSnpIndexConfig = correctAxies.CreateDeltaSnpIndexGraphConfig(_settings);
            var scoreConfig = correctAxies.CreateWindowScoreGraphConfig(_settings);
            var qtlCountConfig = correctAxies.CreateQtlCountGraphConfig(_settings);

            return new GraphCreator[]
            {
                new Bulk1SnpIndexGraphCreator(snpIndexConfig),
                new Bulk2SnpIndexGraphCreator(snpIndexConfig),
                new BulkSnpIndexGraphCreator(snpIndexConfig),
                new DeltaSnpIndexGraphCreator(deltaSnpIndexConfig),
                new WindowScoreGraphCreator(scoreConfig),
                new QtlVariantCountGraphCreator(qtlCountConfig)
            };
        }
    }
}
