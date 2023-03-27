using PolyploidQtlSeqCore.Options;

namespace PolyploidQtlSeqCore.QtlAnalysis.OxyGraph
{
    /// <summary>
    /// X軸の目盛り間隔(MB)
    /// </summary>
    public class XAxisMajorStep
    {
        /// <summary>
        /// 目盛り間隔の最小値
        /// </summary>
        public const int MINIMUM = 1;

        /// <summary>
        /// 目盛り間隔の最大値
        /// </summary>
        public const int MAXIMUM = 50;

        /// <summary>
        /// 目盛り間隔の規定値
        /// </summary>
        public const int DEFAULT = 5;

        /// <summary>
        /// オプションスイッチのShortName
        /// </summary>
        public const string SHORT_NAME = "xs";

        /// <summary>
        /// オプションスイッチのLongName
        /// </summary>
        public const string LONG_NAME = "xStep";

        /// <summary>
        /// オプションスイッチの説明
        /// </summary>
        public const string DESCRIPTION = "X axis major step (MB).";

        /// <summary>
        /// データ検証エラーメッセージ
        /// </summary>
        public const string VALIDATION_ERROR_MESSAGE = "The -xs option must be an integer greater than 1 and less than or equal to 50.";

        /// <summary>
        /// X軸の目盛り間隔(MB)を作成する。
        /// </summary>
        /// <param name="value">ステップ数(MB)</param>
        /// <param name="parameterDictionary">パラメータファイルの中身</param>
        /// <param name="userOptionDictionary">ユーザー指定LongName辞書</param>
        public XAxisMajorStep(int value, IReadOnlyDictionary<string, string> parameterDictionary,
            IReadOnlyDictionary<string, bool> userOptionDictionary)
        {
            Value = OptionValue.GetValue(LONG_NAME, value, parameterDictionary, userOptionDictionary);

            if (Value < MINIMUM || Value > MAXIMUM) throw new ArgumentException(VALIDATION_ERROR_MESSAGE);
        }

        /// <summary>
        /// X軸の目盛り間隔(MB)を取得する。
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
