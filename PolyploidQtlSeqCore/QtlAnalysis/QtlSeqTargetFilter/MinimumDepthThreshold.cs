using PolyploidQtlSeqCore.Options;

namespace PolyploidQtlSeqCore.QtlAnalysis.QtlSeqTargetFilter
{
    /// <summary>
    /// (全サンプルの）最低Depthしきい値
    /// </summary>
    public class MinimumDepthThreshold
    {
        /// <summary>
        /// 最低Depthしきい値の最小値
        /// </summary>
        public const int MINIMUM = 1;

        /// <summary>
        /// 最低Depthしきい値の最大値
        /// </summary>
        public const int MAXIMUM = 10000;

        /// <summary>
        /// 最低Depthしきい値の規定値
        /// </summary>
        public const int DEFAULT = 40;

        /// <summary>
        /// オプションスイッチのShortName
        /// </summary>
        public const string SHORT_NAME = "md";

        /// <summary>
        /// オプションスイッチのLongName
        /// </summary>
        public const string LONG_NAME = "minDepth";

        /// <summary>
        /// オプションスイッチの説明
        /// </summary>
        public const string DESCRIPTION = "Minimum Depth threshold. The variants with even one sample below this threshold are excluded for QTL analysis.";

        /// <summary>
        /// データ検証エラーメッセージ
        /// </summary>
        public const string VALIDATION_ERROR_MESSAGE = "The -md option must be an integer greater than or equal to 1 and less than or equal to 10000.";

        /// <summary>
        /// 最低Depthしきい値を作成する。
        /// </summary>
        /// <param name="threshold">しきい値</param>
        /// <param name="parameterDictionary">パラメータファイルの中身</param>
        /// <param name="userOptionDictionary">ユーザー指定LongName辞書</param>
        public MinimumDepthThreshold(int threshold, IReadOnlyDictionary<string, string> parameterDictionary,
            IReadOnlyDictionary<string, bool> userOptionDictionary)
        {
            Value = OptionValue.GetValue(LONG_NAME, threshold, parameterDictionary, userOptionDictionary);

            if (Value < MINIMUM || Value > MAXIMUM) throw new ArgumentException(VALIDATION_ERROR_MESSAGE);
        }

        /// <summary>
        /// 最低Depthしきい値を取得する。
        /// </summary>
        internal int Value { get; }

        /// <summary>
        /// パラメータファイル記載用行テキストに変換する。
        /// </summary>
        /// <returns>パラメータ行テキスト</returns>
        internal string ToParameterFileLine()
        {
            return $"{LONG_NAME}\t{Value}";
        }

    }
}
