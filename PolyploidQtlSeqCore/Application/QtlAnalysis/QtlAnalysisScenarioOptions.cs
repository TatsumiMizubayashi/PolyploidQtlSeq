using PolyploidQtlSeqCore.QtlAnalysis;
using PolyploidQtlSeqCore.QtlAnalysis.Distribution;
using PolyploidQtlSeqCore.QtlAnalysis.IO;
using PolyploidQtlSeqCore.QtlAnalysis.OxyGraph;
using PolyploidQtlSeqCore.QtlAnalysis.QtlSeqTargetFilter;
using PolyploidQtlSeqCore.QtlAnalysis.SlidingWindow;
using PolyploidQtlSeqCore.Share;

namespace PolyploidQtlSeqCore.Application.QtlAnalysis
{
    /// <summary>
    /// QTL解析シナリオオプション
    /// </summary>
    internal class QtlAnalysisScenarioOptions
    {
        private static readonly IReadOnlyDictionary<string, string> _toLongNameDictionary;

        /// <summary>
        /// コンストラクタ
        /// </summary>
        static QtlAnalysisScenarioOptions()
        {
            var toLongNameDictionary = new Dictionary<string, string>()
            {
                [OutputDirectory.SHORT_NAME] = OutputDirectory.LONG_NAME,
                [OutputDirectory.LONG_NAME] = OutputDirectory.LONG_NAME,

                [DisplayAnnotationImpacts.SHORT_NAME] = DisplayAnnotationImpacts.LONG_NAME,
                [DisplayAnnotationImpacts.LONG_NAME] = DisplayAnnotationImpacts.LONG_NAME,

                [ThreadNumber.SHORT_NAME] = ThreadNumber.LONG_NAME,
                [ThreadNumber.LONG_NAME] = ThreadNumber.LONG_NAME
            };

            QtlSeqTargetPolicyOption.AddLongNameKeyValuePair(toLongNameDictionary);
            NoQtlDistributionOption.AddLongNameKeyValuePair(toLongNameDictionary);
            SlidingWindowAnalysisOption.AddLongNameKeyValuePair(toLongNameDictionary);
            GraphOption.AddLongNameKeyValuePair(toLongNameDictionary);

            _toLongNameDictionary = toLongNameDictionary;
        }

        /// <summary>
        /// LongName変換項目を辞書に追加する。
        /// </summary>
        /// <param name="dictionary">LongName変換辞書</param>
        public static void AddLongNameKeyValuePair(Dictionary<string, string> dictionary)
        {
            foreach (var keyValuePair in _toLongNameDictionary)
            {
                dictionary.Add(keyValuePair.Key, keyValuePair.Value);
            }
        }

        /// <summary>
        /// QTL解析シナリオオプションを作成する。
        /// </summary>
        /// <param name="optionValues">オプションの値</param>
        /// <param name="parameterDictionary">パラメータファイルの中身</param>
        /// <param name="userOptionDictionary">ユーザー指定オプション辞書</param>
        public QtlAnalysisScenarioOptions(IQtlAnalysisScenarioOptions optionValues, IReadOnlyDictionary<string, string> parameterDictionary,
            IReadOnlyDictionary<string, bool> userOptionDictionary)
        {
            OutputDir = new OutputDirectory(optionValues.OutputDir, parameterDictionary, userOptionDictionary);
            DisplayAnnotationImpacts = new DisplayAnnotationImpacts(optionValues.DisplayAnnotationImpacts, parameterDictionary, userOptionDictionary);
            ThreadNumber = new ThreadNumber(optionValues.ThreadNumber, parameterDictionary, userOptionDictionary);

            QtlSeqTargetPolicyOption = new QtlSeqTargetPolicyOption(optionValues, parameterDictionary, userOptionDictionary);
            NoQtlDistributionOption = new NoQtlDistributionOption(optionValues, parameterDictionary, userOptionDictionary);
            SlidingWindowAnalysisOption = new SlidingWindowAnalysisOption(optionValues, parameterDictionary, userOptionDictionary);
            GraphOption = new GraphOption(optionValues, parameterDictionary, userOptionDictionary);
        }

        /// <summary>
        /// 出力ディレクトリ
        /// </summary>
        public OutputDirectory OutputDir { get; }

        /// <summary>
        /// 表示するAnnotation Impact
        /// </summary>
        public DisplayAnnotationImpacts DisplayAnnotationImpacts { get; }

        /// <summary>
        /// スレッド数
        /// </summary>
        public ThreadNumber ThreadNumber { get; }

        /// <summary>
        /// QTL-seq対象ポリシーオプション
        /// </summary>
        public QtlSeqTargetPolicyOption QtlSeqTargetPolicyOption { get; }

        /// <summary>
        /// QTLなし分布オプション
        /// </summary>
        public NoQtlDistributionOption NoQtlDistributionOption { get; }

        /// <summary>
        /// SlidingWindow解析オプション
        /// </summary>
        public SlidingWindowAnalysisOption SlidingWindowAnalysisOption { get; }

        /// <summary>
        /// グラフオプション
        /// </summary>
        public GraphOption GraphOption { get; }

        /// <summary>
        /// パラメータファイルに記載する行テキストに変換する。
        /// </summary>
        /// <returns>パラメータ行テキスト</returns>
        public string[] ToParameterFileLines()
        {
            return new[] { OutputDir.ToParameterFileLine() }
                .Concat(QtlSeqTargetPolicyOption.ToParameterFileLines())
                .Concat(NoQtlDistributionOption.ToParameterFileLines())
                .Concat(SlidingWindowAnalysisOption.ToParameterFileLines())
                .Concat(GraphOption.ToParameterFileLines())
                .Concat(new[]
                {
                    DisplayAnnotationImpacts.ToParameterFileLine(),
                    ThreadNumber.ToParameterFileLine()
                })
                .ToArray();
        }
    }
}
