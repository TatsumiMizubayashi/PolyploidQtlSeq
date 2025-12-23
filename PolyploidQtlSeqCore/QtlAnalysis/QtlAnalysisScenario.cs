using Kurukuru;
using PolyploidQtlSeqCore.QtlAnalysis.OxyGraph;
using PolyploidQtlSeqCore.QtlAnalysis.Preprocess;
using PolyploidQtlSeqCore.QtlAnalysis.IO;
using PolyploidQtlSeqCore.QtlAnalysis.Preprocess.IO;
using PolyploidQtlSeqCore.QtlAnalysis.SlidingWindow;

namespace PolyploidQtlSeqCore.QtlAnalysis
{
    /// <summary>
    /// QTL解析シナリオ
    /// </summary>
    internal class QtlAnalysisScenario
    {
        private readonly QtlAnalysisScenarioSettings _settings;

        /// <summary>
        /// QTL解析シナリオを作成する。
        /// </summary>
        /// <param name="settings">設定</param>
        public QtlAnalysisScenario(QtlAnalysisScenarioSettings settings)
        {
            _settings = settings;
        }

        /// <summary>
        /// QTL解析を実行する。
        /// </summary>
        /// <param name="inputVcf">入力VCFファイル</param>
        /// <returns>終了コード</returns>
        public int Run(InputVcf inputVcf)
        {
            var code = 0;
            Spinner.Start("Loading VCF file ...", spinner =>
            {
                try
                {
                    _settings.OutputDir.Create();

                    var snpIndexVariants = GetTargetSnpIndexVariants(inputVcf);

                    spinner.Text = "Variant QTL analysis ...";
                    var qtlVariants = AnalyzeVariantQtl(snpIndexVariants);

                    spinner.Text = "Sliding window QTL analysis ...";
                    var windows = AnalyzeSlidingWindow(qtlVariants);
                    var swQtlVariants = LinkSlidingWindowQtl(qtlVariants, windows);

                    spinner.Text = "Saving ...";
                    SnpIndexFileCreator.Create(_settings.OutputDir, swQtlVariants, _settings.DisplayAnnotationImpacts);
                    SlidingWindowFileCreator.Create(_settings.OutputDir, windows);
                    var graphCreator = new QtlSeqGraphCreator(_settings.GraphSettings);
                    graphCreator.Create(_settings.OutputDir, swQtlVariants, windows);

                    spinner.Succeed("QTL-Seq analysis completed");
                }
                catch (Exception ex)
                {
                    code = 1;
                    spinner.Fail("QTL-Seq analysis error");
                    Console.Error.WriteLine(ex.Message);
                    Console.Error.WriteLine(ex.StackTrace);
                }
            });

            return code;
        }

        /// <summary>
        /// VCFファイルを読み込み、QTL解析対象変異情報を取得する。
        /// </summary>
        /// <param name="inputVcf">入力VCFファイル</param>
        /// <returns>解析対象変異情報</returns>
        private SnpIndexVariant[] GetTargetSnpIndexVariants(InputVcf inputVcf)
        {
            var analyzableVriantPolicy = new AnalyzableVariantPolicy();
            var qtlSeqTargetVariantPolicy = _settings.QtlSeqTargetPolicySettings.CreatePolicy();

            var vcfVariants = VcfFileParser.Parse(inputVcf.Path);
            var variants = vcfVariants
                .Where(x => analyzableVriantPolicy.ComplyWithAll(x))
                .Select(x => new SnpIndexVariant(x))
                .Where(x => qtlSeqTargetVariantPolicy.ComplyWithAll(x))
                .ToArray();
            if (variants.Length == 0) 
                throw new InvalidOperationException("No variants are available for analysis. Try loosening the filter criteria. If that still doesn’t resolve the issue, verify that the VCF file contains variants that meet the required analysis conditions.");

            return variants;
        }

        /// <summary>
        /// 変異のQTL解析を行う。
        /// </summary>
        /// <param name="variants">解析対象変異</param>
        /// <returns>QTL解析済み変異</returns>
        private SnpIndexVariantWithQtl[] AnalyzeVariantQtl(SnpIndexVariant[] variants)
        {
            var analyzer = new VariantQtlAanlyzer(_settings.NoQtlDistributionSettings, _settings.ThreadNumber);
            return analyzer.Analyze(variants);
        }

        /// <summary>
        /// Sliding Window解析を行う。
        /// </summary>
        /// <param name="variants">変異情報</param>
        /// <returns>SlidingWindow</returns>
        private Window[] AnalyzeSlidingWindow(SnpIndexVariantWithQtl[] variants)
        {
            var slidingWindowAnalyzer = new SlidingWindowAnalyzer(_settings.SlidingWindowAnalysisSettings, _settings.ThreadNumber);
            return slidingWindowAnalyzer.Analyze(variants);
        }

        /// <summary>
        /// 変異情報にSlidingWindow解析結果を紐づける。
        /// </summary>
        /// <param name="variants">変異情報</param>
        /// <param name="windows">SlidingWindow</param>
        /// <returns>SlidingWindow結果紐づけ済み変異情報</returns>
        private SnpIndexVariantWithSlidingWindowQtl[] LinkSlidingWindowQtl(SnpIndexVariantWithQtl[] variants, Window[] windows)
        {
            var creator = new SnpIndexVariantWithSlidingWindowQtlCreator(windows, _settings.ThreadNumber);
            return creator.Create(variants);
        }

    }
}
