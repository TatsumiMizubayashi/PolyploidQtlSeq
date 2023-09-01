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
        private readonly QtlAnalysisScenarioOptions _option;

        /// <summary>
        /// QTL解析シナリオを作成する。
        /// </summary>
        /// <param name="option">オプション</param>
        public QtlAnalysisScenario(QtlAnalysisScenarioOptions option)
        {
            _option = option;
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
                    _option.OutputDir.Create();

                    var snpIndexVariants = GetTargetSnpIndexVariants(inputVcf);

                    spinner.Text = "Variant QTL analysis ...";
                    var qtlVariants = AnalyzeVariantQtl(snpIndexVariants);

                    spinner.Text = "Sliding window QTL analysis ...";
                    var windows = AnalyzeSlidingWindow(qtlVariants);
                    var swQtlVariants = LinkSlidingWindowQtl(qtlVariants, windows);

                    spinner.Text = "Saving ...";
                    SnpIndexFileCreator.Create(_option.OutputDir, swQtlVariants, _option.DisplayAnnotationImpacts);
                    SlidingWindowFileCreator.Create(_option.OutputDir, windows);
                    var graphCreator = new QtlSeqGraphCreator(_option.GraphOption);
                    graphCreator.Create(_option.OutputDir, swQtlVariants, windows);

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
            var qtlSeqTargetVariantPolicy = _option.QtlSeqTargetPolicyOption.CreatePolicy();

            var vcfVariants = VcfFileParser.Parse(inputVcf.Path);
            return vcfVariants
                .Where(x => analyzableVriantPolicy.ComplyWithAll(x))
                .Select(x => new SnpIndexVariant(x))
                .Where(x => qtlSeqTargetVariantPolicy.ComplyWithAll(x))
                .ToArray();
        }

        /// <summary>
        /// 変異のQTL解析を行う。
        /// </summary>
        /// <param name="variants">解析対象変異</param>
        /// <returns>QTL解析済み変異</returns>
        private SnpIndexVariantWithQtl[] AnalyzeVariantQtl(SnpIndexVariant[] variants)
        {
            var analyzer = new VariantQtlAanlyzer(_option.NoQtlDistributionOption, _option.ThreadNumber);
            return analyzer.Analyze(variants);
        }

        /// <summary>
        /// Sliding Window解析を行う。
        /// </summary>
        /// <param name="variants">変異情報</param>
        /// <returns>SlidingWindow</returns>
        private Window[] AnalyzeSlidingWindow(SnpIndexVariantWithQtl[] variants)
        {
            var slidingWindowAnalyzer = new SlidingWindowAnalyzer(_option.SlidingWindowAnalysisOption, _option.ThreadNumber);
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
            var creator = new SnpIndexVariantWithSlidingWindowQtlCreator(windows, _option.ThreadNumber);
            return creator.Create(variants);
        }

    }
}
