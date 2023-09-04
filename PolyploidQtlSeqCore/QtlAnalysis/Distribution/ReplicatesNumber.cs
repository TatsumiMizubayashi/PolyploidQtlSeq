using PolyploidQtlSeqCore.Options;

namespace PolyploidQtlSeqCore.QtlAnalysis.Distribution
{
    /// <summary>
    /// 分布作成時の試行回数
    /// </summary>
    public class ReplicatesNumber
    {
        /// <summary>
        /// 試行回数の最小値
        /// </summary>
        private const int MINIMUM = 1000;

        /// <summary>
        /// 試行回数の最大値
        /// </summary>
        private const int MAXIMUM = 1_000_000;

        /// <summary>
        /// 試行回数の規定値
        /// </summary>
        [Obsolete("削除予定")]
        public const int DEFAULT = 5000;

        /// <summary>
        /// オプションスイッチのShortName
        /// </summary>
        [Obsolete("削除予定")]
        public const string SHORT_NAME = "N";

        /// <summary>
        /// オプションスイッチのLongName
        /// </summary>
        [Obsolete("削除予定")]
        public const string LONG_NAME = "NRep";

        /// <summary>
        /// オプションスイッチの説明
        /// </summary>
        [Obsolete("削除予定")]
        public const string DESCRIPTION = "Number of simulation replicates to generate a null distribution which is free from QTLs.";

        /// <summary>
        /// データ検証エラーメッセージ
        /// </summary>
        [Obsolete("削除予定")]
        public const string VALIDATION_ERROR_MESSAGE = "The -N option must be an integer greater than or equal to 1000 and less than or equal to 1,000,000.";

        /// <summary>
        /// 分布作成時の試行回数を作成する。
        /// </summary>
        /// <param name="number">試行回数</param>
        /// <param name="parameterDictionary">パラメータファイルの中身</param>
        /// <param name="userOptionDictionary">ユーザー指定LongName辞書</param>
        public ReplicatesNumber(int number, IReadOnlyDictionary<string, string> parameterDictionary,
            IReadOnlyDictionary<string, bool> userOptionDictionary)
        {
            Value = OptionValue.GetValue(LONG_NAME, number, parameterDictionary, userOptionDictionary);

            if (Value < MINIMUM || Value > MAXIMUM) throw new ArgumentException(VALIDATION_ERROR_MESSAGE);
        }

        /// <summary>
        /// 分布作成時の試行回数を作成する。
        /// </summary>
        /// <param name="number">試行回数</param>
        public ReplicatesNumber(int number)
        {
            if (number < MINIMUM || number > MAXIMUM) throw new ArgumentOutOfRangeException(nameof(number));

            Value = number;
        }


        /// <summary>
        /// Bulk2個体数を取得する。
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
